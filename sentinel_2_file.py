#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import snappy
import numpy
from snappy import ProductIO
import matplotlib.image as mpimg
import os
from mpl_toolkits.basemap import Basemap
import argparse
import sys
import glob
from snappy import GPF, Mask
import json
import rasterio.features
import shapely.geometry
from PIL import Image
import polymer
import polymer.level2
import polymer.ancillary
import polymer.level1
import polymer.main
import cPickle as pickle
import netCDF4 as nc


jpy = snappy.jpy

masks = ['lc_cloud', 'lc_cloud_sure', 'lc_cloud_buffer', 'lc_cloud_shadow', 'lc_cirrus_sure', 'lc_cirrus_ambiguous', 'lc_clear_land', 'lc_clear_water']

class sentinel2_msi_file(object):
    def __init__(self,sentinel_2_file_path, chunk_size=2000):
        self.full_path = sentinel_2_file_path
        self.source_filename = os.path.basename(sentinel_2_file_path)
        print self.source_filename
        self.product = ProductIO.readProduct(sentinel_2_file_path)
        self.atm_product = None
        self.atm_band_names=None
        self.width = self.product.getSceneRasterWidth()
        self.height = self.product.getSceneRasterHeight()
        self.name = self.product.getName()
        self.description = self.product.getDescription()
        self.band_names = self.product.getBandNames()
        self.chunk_size = chunk_size
        self.chunks = self.define_tiles(self.width, self.height)
        self.ten_m_bands = ['B2', 'B3', 'B4', 'B8']
        self.twenty_m_bands = ['B5', 'B6', 'B7', 'B8A', 'B11', 'B12']
        self.sixty_m_bands = ['B1', 'B9', 'B10']
        self.tile = None
        self.true_chunk_def = [0,0]
        self.current_tile = None
        self.tile_north = None
        self.tile_south = None
        self.tile_east = None
        self.tile_west = None
        self.tile_water = False
        print "doing atm corr"
        self.get_atm_corr_file()
        print "reading as rgb image"
        self.make_rgb_data()
        self.resamp_geocode = None
        self.do_resample()
        self.do_idepix()
        self.json = {"labels" : [], "complete": None}
        print "done"
    
    def get_correct_sizes(self, band_name, x, y):
        if band_name in self.ten_m_bands or band_name == None:
            chunk_size = self.chunk_size
            x = x
            y = y
        elif band_name in self.twenty_m_bands:
            chunk_size = self.chunk_size / 2
            x = x / 2
            y = y / 2
        elif band_name in self.sixty_m_bands:
            chunk_size = self.chunk_size / 6
            x = x / 6
            y = y / 6
        else:
            print band_name
        return chunk_size, x, y

    def read_tile(self, x, y):
        bands = {}
        band_names = [band for band in self.band_names if band.startswith('B')]
        for band_name in band_names:
            chunk_size, startx, starty = self.get_correct_sizes(None, x, y)
            band = self.resamp.getBand(band_name)
            band_array = numpy.zeros((chunk_size, chunk_size), dtype=numpy.float32)
            band.readPixels(startx,starty,chunk_size,chunk_size,band_array)
            bands[band_name] = band_array
        
        atm_band_names = [band for band in self.atm_band_names if band.startswith('R') or band.startswith('log')]
        for band_name in atm_band_names:
            chunk_size, startx, starty = self.get_correct_sizes(None, x, y)
            band = self.atm_resampled_product.getBand(band_name)
            band_array = numpy.zeros((chunk_size, chunk_size), dtype=numpy.float32)
            band.readPixels(startx,starty,chunk_size,chunk_size,band_array)
            bands[band_name] = band_array
        return bands

    def define_tiles(self, width, height):
        x = [(x-self.chunk_size, x-1) for x in range(self.chunk_size, width, self.chunk_size)]
        y = [(y-self.chunk_size, y-1) for y in range(self.chunk_size, height, self.chunk_size)]

        if x[-1][1] != width:
            x.append((x[-1][1] + 1, width))
        
        if y[-1][1] != height:
            y.append((y[-1][1] + 1, height))
        chunks = []
        for w in x:
            for h in y:
                chunks.append([w,h])
        return chunks

    def get_tile(self, tile_no):
        chunk = self.chunks[tile_no]
        tile = self.read_tile(chunk[0][0], chunk[1][0])
        self.tile = tile
        self.current_tile = tile_no
        gc = self.product.getSceneGeoCoding()
        print gc
        bot_pos = gc.getGeoPos(snappy.PixelPos(chunk[0][1], chunk[1][1]),None)
        top_pos = gc.getGeoPos(snappy.PixelPos(chunk[0][0], chunk[1][0]),None)
        self.tile_north = top_pos.getLat()
        self.tile_east = bot_pos.getLon()
        self.tile_south = bot_pos.getLat()
        self.tile_west = top_pos.getLon()
        self.lats = numpy.arange(self.tile_south, self.tile_north, ((self.tile_north - self.tile_south) / (self.chunk_size-1)))
        self.lons = numpy.arange(self.tile_east, self.tile_west, ((self.tile_east - self.tile_west) / (self.chunk_size-1)))
        self.test_tile_water()


    def test_tile_water(self):
        self.tile_water = False
        try:
           self.basemap = Basemap(llcrnrlon=self.tile_west,
                                  llcrnrlat=self.tile_south,
                                  urcrnrlon=self.tile_east,
                                  urcrnrlat=self.tile_north,
                                  resolution='h',
                                  projection='tmerc',
                                  lon_0=self.tile_west + self.tile_east / 2,
                                  lat_0=self.tile_south + self.tile_north / 2)
           self.basemap.drawlsmask(resolution='f', grid=1.25)
           self.tile_water = numpy.any(self.basemap.lsmask == 0)
        except:
           self.tile_water = False
        
    def generate_rgb_tile(self, output_location, r='B4', g='B3', b='B2'):
        img = numpy.dstack((self.tile[r], self.tile[g], self.tile[b]))
        img_max = img.max()
        print img_max
        img = numpy.multiply(img, 255.0)
        img = numpy.divide(img, img_max)
        img = img.astype(numpy.uint8)
        mpimg.imsave(output_location, img, format='jpg')

    def make_rgb_data(self, r='B4', g='B3', b='B2'):
        b = self.product.getBand(b)
        g = self.product.getBand(g)
        r = self.product.getBand(r)
        print "read bands"
        info = snappy.ProductUtils.createImageInfo([r,g,b], True,snappy.ProgressMonitor.NULL)
        print "set info"
        self.rgb_image = snappy.ProductUtils.createRgbImage([r,g,b],info, snappy.ProgressMonitor.NULL)
        print "rgb set"

    def generate_snappy_rgb_tile(self, output_location):
        chunk = self.chunks[self.current_tile]
        print self.rgb_image
        try:
            tile = self.rgb_image.getSubimage(chunk[0][0], chunk[1][0], self.chunk_size, self.chunk_size)
            self.true_chunk_def = [self.chunk_size, self.chunk_size]
        except:
            if chunk[0][0] + self.chunk_size >= self.width:
                self.true_chunk_def[0] = self.width - chunk[0][0]
            if chunk[1][0] + self.chunk_size >= self.height:
                self.true_chunk_def[1] = self.height - chunk[0][0]
            try:
                tile = self.rgb_image.getSubimage(chunk[0][0], chunk[1][0], self.true_chunk_def[0], self.true_chunk_def[1])
            except:
                print chunk[0][0], chunk[1][0], self.true_chunk_def[0], self.true_chunk_def[1]
                print self.width, self.height
                #can't do any thing more
                return True
        RenderedImage = jpy.get_type('java.awt.image.RenderedImage')
        rendered_rgb_image = jpy.cast(tile, RenderedImage)
        File = jpy.get_type('java.io.File')
        oname = File(output_location)
        imageio = jpy.get_type('javax.imageio.ImageIO')
        imageio.write(rendered_rgb_image, "png", oname)
        self.convert_rgb_to_jpg(output_location)
        self.make_tile_mask_json(output_location)
        self.make_text_definitions(output_location)
        self.make_spectral_nc(output_location)

    def make_spectral_nc(self, output_location):
        print "dumping spectral data"
        out_file = output_location.replace("png", "nc")
        spec_data = self.tile
        shape = spec_data[spec_data.keys()[0]].shape
        print out_file
        out_file=nc.Dataset(out_file,"w", format="NETCDF4")

        out_file.createDimension('lon', shape[0])
        out_file.createDimension('lat', shape[1])
        longitude = out_file.createVariable('Longitude', 'f4', 'lon', zlib=True, complevel=9)
        latitude = out_file.createVariable('Latitude', 'f4', 'lat', zlib=True, complevel=9)

        longitude[:] = range(0, shape[0])
        latitude[:] = range(0, shape[1])

        for f in spec_data.keys():
            print f
            try:
                temp = out_file.createVariable(f, 'f4', ('lon', 'lat'), zlib=True, complevel=9, least_significant_digit=3)
            except Exception as e:
                print e
                temp = out_file.variables[f]
            temp[:] = spec_data[f]

        out_file.close()


    def convert_rgb_to_jpg(self, png_file):
        im = Image.open(png_file)
        rgb_im = im.convert('RGB')
        rgb_im.save(png_file.replace(".png", ".jpg"))

    def do_resample(self):
        HashMap = snappy.jpy.get_type('java.util.HashMap')
        snappy.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
        self.parameters = HashMap()
        self.parameters.put('paramName',False)
        self.parameters.put('targetResolution',10)
        result = snappy.GPF.createProduct('Resample', self.parameters, self.product)
        self.resamp_geocode = result.getBand("B3").getGeoCoding()
        self.resamp = result
        
    def do_idepix(self):
        self.idepix_prod = GPF.createProduct('Idepix.Sentinel2',self.parameters,self.resamp)
    
    def export_mask_as_polygon(self, mask_name):
        mask = self.idepix_prod.getMaskGroup().get(mask_name)
        real_mask = jpy.cast(mask, Mask)
        mask_arr = numpy.zeros((self.chunk_size, self.chunk_size))
        chunk = self.chunks[self.current_tile]
        real_mask.readPixels(chunk[0][0], chunk[1][0], self.true_chunk_def[0], self.true_chunk_def[1], mask_arr)
        mask_one_arr = mask_arr
        mask_one_arr[mask_one_arr > 1] = 1
        mask_one_arr = mask_one_arr.astype(numpy.int16)
        shapes = rasterio.features.shapes(mask_one_arr)
        polygons = [shapely.geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] == 1]
        return polygons
    
    def make_text_definitions(self, output_location):
        chunk = self.chunks[self.current_tile]
        nwpp = snappy.PixelPos(chunk[0][0], chunk[1][0])
        sepp = snappy.PixelPos(chunk[0][0] +self.true_chunk_def[0], chunk[1][0] +self.true_chunk_def[1])
        nwgp = self.resamp_geocode.getGeoPos(nwpp, None)
        segp = self.resamp_geocode.getGeoPos(sepp, None)
        bounds = {}
        bounds["pixel"] = {}
        bounds["pixel"]["north"] = chunk[0][0]
        bounds["pixel"]["south"] = chunk[0][0] + self.true_chunk_def[1]
        bounds["pixel"]["east"] = chunk[1][0] + self.true_chunk_def[0]
        bounds["pixel"]["west"] = chunk[1][0]
        bounds["latlon"] = {}
        bounds["latlon"]["north"] = nwgp.getLat()
        bounds["latlon"]["south"] = segp.getLat() 
        bounds["latlon"]["east"] = segp.getLon()
        bounds["latlon"]["west"] = nwgp.getLon()
        print "outputting bounds"
        with open(os.path.join(output_location.replace(".png", "_bounds.json")), "w") as output:
            json.dump(bounds, output)
        return True


    def make_tile_mask_json(self, output_location):
        labels = []
        object_id = 1
        json_output = self.json
        if os.path.isfile(os.path.join(output_location, self.tile_name + "__labels.json")):
            old_json = open(os.path.join(output_location, self.tile_name + "__labels.json")).read()
            old_json = json.loads(old_json)
            json_output = old_json
            labels = old_json['labels']
            if labels:
                for lab in labels:
                    if lab["object_id"] > object_id:
                        object_id = lab["object_id"]
            object_id += 1
            
        for mask in masks:
            polys = self.export_mask_as_polygon(mask)
            for poly in polys:
                coords = []
                for coord in numpy.array(poly.exterior):
                    pixels = {"x" : coord[0], "y" : coord[1]}
                    coords.append(pixels)
                mask_name = self.get_mask_realname(mask)
                labels.append({"label_type" : "polygon", "label_class": mask_name, "vertices" : coords, "object_id" : object_id})
                object_id += 1

        json_output['labels'] = labels
        with open(os.path.join(output_location.replace(".png", "__labels_new.json")), "w") as output:
            json.dump(json_output, output)
        

    def get_mask_realname(self, mask):
        if "cloud" in mask:
            return "cloud"
        if "water" in mask:
            return "water"
        if "land" in mask:
            return "land"

    def get_atm_corr_file(self, atm_corr_storage="/data/gfs03/scratch/stgo/atm_corrected_sen2", ancillary_dir="/data/datasets/operational/ancillary/obpg/live_mirror/Meteorological/"):
        print self.source_filename
        print self.full_path
        atm_file = os.path.join(atm_corr_storage, self.source_filename + "_atm.nc")
        print atm_file
        if not os.path.isfile(atm_file):
            self.ancillary = polymer.ancillary.Ancillary_NASA(directory=ancillary_dir, 
                                                         offline=True)
            self.lev1 = polymer.level1.Level1(self.full_path,
                                         sensor="msi",
                                         ancillary=self.ancillary)
            self.lev2 = polymer.level2.Level2(filename=atm_file,
                                         fmt="netcdf4",
                                         datasets=polymer.level2.default_datasets.extend(polymer.level2.analysis_datasets))

            polymer.main.run_atm_corr(self.lev1,
                        self.lev2,
                        initial_point_1=[0,0],
                        BITMASK_INVALID=4+512,
                        min_abs=0,
                        initial_point_2=[1, 1],
                        thres_Rcloud=-1,
                        thres_Rcloud_std=-1,
                        reinit_rw_neg= 1,
                        water_model="PR05",
                        BITMASK_REJECT=0,
                        bands_rw=[443,490,560,665,705,740,783,842,865,1375,1610,2190]
                        )

        atm_product = ProductIO.readProduct(atm_file)
        HashMap = snappy.jpy.get_type('java.util.HashMap')
        snappy.GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
        parameters = HashMap()
        parameters.put('targetWidth', self.width)
        parameters.put('targetHeight', self.height)
        self.atm_resampled_product = GPF.createProduct('Resample', parameters, atm_product)
        self.atm_band_names = self.atm_resampled_product.getBandNames()
                   

        return atm_file

    def create_rgb_for_all_tiles(self, output_location, r='B4', g='B3', b='B2', test_water=False):
        out_file_name = self.source_filename + "_tile_{}.png"
        for tile in range(0, len(self.chunks)):
            print "creating tile {}".format(tile)
            tile_name = out_file_name.format(tile)
            self.tile_name = tile_name
            self.json["image_filename"] = self.tile_name
            self.get_tile(tile)
            print self.tile_water
            if self.tile_water or not test_water:
                print "making", tile_name
                self.generate_snappy_rgb_tile(os.path.join(output_location, tile_name))

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description='create some sentinel jpgs')
   parser.add_argument('output', help="output location")
   parser.add_argument('filename', help="sentinel 2 folder")
   args = parser.parse_args()
   sen_file = sentinel2_msi_file(args.filename) 
   sen_file.create_rgb_for_all_tiles(args.output)
