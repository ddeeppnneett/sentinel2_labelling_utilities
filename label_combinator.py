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
import shapely.ops
from PIL import Image
import time
import sys
import threading
import multiprocessing
from itertools import product
import cPickle as pickle
import netCDF4 as nc

"""
At some point this should be extended further to push all the json/polygons and rgb data into a single netcdf for distribution purposes.
"""

def fix_polygons(label_file):
    """
    Combines manual and automatically generated labels into one big file, with simplified overlapping polygons.
    """
    try:
        with open(label_file) as jfile:
            base_json = json.load(jfile)
        polygons = base_json['labels']
    except:
        polygons = []
        base_json = None
    if os.path.exists(label_file.replace("__labels.json", "__labels_new.json")):
        with open(label_file.replace("__labels.json", "__labels_new.json")) as jfile:
            second_base_json = json.load(jfile)
        polygons.extend(second_base_json['labels'])
        if not base_json:
            base_json = second_base_json
        
    if not base_json:
        sys.exit(0)
    skip = []
    final_polygons = []
    label_names = set([x['label_class'] for x in polygons])
    objectid = 1
    for label in label_names:
        if label == None:
            continue
        label_polygons = [x for x in polygons if label == x['label_class'] and x['label_type'] == 'polygon']
        label_boxes = [x for x in polygons if label == x['label_class'] and x['label_type'] == 'box']
        final_boxes = []
        for box in label_boxes:
            if box['size']['y'] > 1970 and box['size']['x'] > 1970:
                box_dict = {'vertices' : [{'x' : 0, 'y' : 0}, {'x' : 2000, 'y' : 0}, {'x' : 2000, 'y' : 2000}, {'x' : 2000, 'y' : 0}]}
                final_boxes.append(box_dict)
            else:
                box_dict = {'vertices' : [{'x' : box['centre']['x'] - (box['size']['x'] / 2), 'y' : box['centre']['y'] - (box['size']['y'] / 2)},
                                          {'x' : box['centre']['x'] + (box['size']['x'] / 2), 'y' : box['centre']['y'] - (box['size']['y'] / 2)},
                                          {'x' : box['centre']['x'] + (box['size']['x'] / 2), 'y' : box['centre']['y'] + (box['size']['y'] / 2)},
                                          {'x' : box['centre']['x'] - (box['size']['x'] / 2), 'y' : box['centre']['y'] + (box['size']['y'] / 2)}]}
                final_boxes.append(box_dict)
        label_polygons.extend(final_boxes)
        polys = []
        for base_index, base_poly_dict in enumerate(label_polygons):
            base_poly = [shapely.geometry.Point(x['x'],x['y']) for x in base_poly_dict['vertices']]
            base_poly = shapely.geometry.MultiPoint(base_poly).convex_hull
            polys.append(base_poly)
        union = shapely.ops.cascaded_union(polys)
        if label == None:
            print "*******************************"
        if type(union) == shapely.geometry.multipolygon.MultiPolygon:
            for polygon in union:
                new_poly_dict = {"label_type": "polygon",
                            "label_class": label, 
                            "vertices": [{"x" : coord[0], "y" : coord[1]} for coord in numpy.array(polygon.exterior)], 
                            "object_id": objectid}
                final_polygons.append(new_poly_dict)
                objectid += 1
        else:
            polygon = union
            new_poly_dict = {"label_type": "polygon",
                        "label_class": label, 
                        "vertices": [{"x" : coord[0], "y" : coord[1]} for coord in numpy.array(polygon.exterior)], 
                        "object_id": objectid}
            final_polygons.append(new_poly_dict)
            objectid += 1
        """
        for base_index, base_poly_dict in enumerate(label_polygons):
            print base_index
            base_poly = [shapely.geometry.Point(x['x'],x['y']) for x in base_poly_dict['vertices']]
            base_poly = shapely.geometry.MultiPoint(base_poly).convex_hull
            join = [base_poly]
            if base_index in skip:
                continue
            for comp_index, comp_poly_dict in enumerate(polygons):
                if comp_index == base_index or comp_index in skip:
                    #skip
                    continue
                comp_poly =  [shapely.geometry.Point(x['x'],x['y']) for x in comp_poly_dict['vertices']]
                comp_poly = shapely.geometry.MultiPoint(comp_poly).convex_hull
                if comp_poly.area < 20:
                    continue
                if comp_poly.intersects(base_poly):
                    join.append(comp_poly)
                    skip.append(comp_index)
            try:
                union = shapely.ops.cascaded_union(join)
                base_poly_dict['vertices'] = [{"x" : coord[0], "y" : coord[1]} for coord in numpy.array(union.exterior)]
            except Exception as e:
                continue
            
            if union.area < 20:
                continue

            final_polygons.append(base_poly_dict)
        """
    for index in range(0, len(final_polygons)):
        final_polygons[index]['object_id'] = index

    base_json['labels'] = final_polygons
    with open(os.path.join(label_file.replace("__labels.json", "__labels_final.json")), "w") as output:
        json.dump(base_json, output)
    



def dimensioniser(low, high):
    """
    legacy
    """
    e = high
    w = low
    perc = (((e - w) / 400) % 1)
    change = (400 - (400 * perc)) / 2
    east_change = change
    west_change = change
    dim_error = False
    try:
        too_big = False
        if e + east_change > 2000:
            west_change = west_change + (east_change - (2000 - e))
            east_change = 2000 - e
            too_big =True
        else:
            east_change = change
        if w - change < 0:
            if too_big:
                raise Exception
            west_change = w
            east_change = change + (change - west_change)
            if e + east_change > 2000:
                raise Exception
        else:
            west_change = change
        east = (e + east_change)
        west = (w - west_change)
        no_of_tiles = (east - west) / 400
    except:
        west = 0
        east = 2000
        dim_error = True
        no_of_tiles = (2000 / 400) + 1
    return  west, east, dim_error, no_of_tiles

def tileifier(polygon,poly_type):
    """
    legacy tester of data intersection.
    """
    w,n,e,s= polygon.bounds
    north, south, dim_error_ns, no_of_tiles_ns = dimensioniser(n, s)
    west, east, dim_error_we, no_of_tiles_we = dimensioniser(w, e)
    total_tiles = no_of_tiles_ns * no_of_tiles_we
    tiles = []
    for updown in range(1, int(no_of_tiles_ns)):
        for leftright in range(1, int(no_of_tiles_we)):
            if not leftright == no_of_tiles_we:
                l = ((west - 400) + (400 * leftright))
                r = (west + (400 * leftright))
            elif not dim_error_we:
                l = ((west - 400) + (400 * leftright))
                r = (west + (400 * leftright))
            else:
                l = (east - 400)
                r = (east)
            if not updown == no_of_tiles_ns:
                u = ((north - 400) + (400 * updown))
                d = (north + (400 * updown))
            elif not dim_error_ns:
                u = ((north - 400) + (400 * updown))
                d = (north + (400 * updown))   
            else:
                u = (south - 400)
                d = (south)
            tile = [(u, l), (u, r), (d,r), (d,l)]
            true_tile = [l, u, r, d]
            tile_poly = shapely.geometry.MultiPoint([shapely.geometry.Point(x) for x in tile]).convex_hull
            intersect = (tile_poly.intersection(polygon).area/tile_poly.area)*100
            inter_min = 0
            inter_max = 101

            if poly_type == 'land':
                inter_max = 80

            if poly_type == 'water':
                inter_min = 70
            
            if poly_type == 'algal':
                inter_min = 80
            
            if poly_type == 'cloud':
                inter_min = 70

            if intersect > inter_min and intersect < inter_max:
                tiles.append((true_tile, tile_poly.intersection(polygon)))
    return tiles

def make_spectral_nc(spec_data, oname):
    """
    Copied from sentinel_2_file and muddled a bit. Outputs a netcdf of a dict.
    """
    print "dumping spectral data"
    shape = spec_data[spec_data.keys()[0]].shape
    print oname
    out_file=nc.Dataset(oname,"w", format="NETCDF4")

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

def image_slicer(arguments):
    """
    Over built driver function, slices tiles up even further so they can be handled by low memory system.
    """
    image_file, output_location = arguments
    spectral_data = None
    print image_file.replace("png","nc")
    if os.path.exists(image_file.replace("png","nc")):
        try:
            nc_file = image_file.replace("png","nc")
            print "trying to read"
            ncd = nc.Dataset(nc_file)
            print nc_file
            print ncd
            spectral_data = {ncd.variables[x].name: ncd.variables[x][:] for x in ncd.variables if len(ncd.variables[x].dimensions) == 2}
        except Exception as e:
            print e
            spectral_data = None

    fix_polygons(image_file.replace(".png", "__labels.json"))
    image_json = image_file.replace(".png", "__labels_final.json")
    with open(image_json) as jfile:
        labels_json = json.load(jfile)
    image = Image.open(image_file)
    labels = [x for x in labels_json['labels']]
    super_json = []
    image_tiles = [[j,i,j+400, i+400] for i in range(0,2000,400) for j in range(0,2000,400)]
    segmentations = []
    for index, tile in enumerate(image_tiles):
        skip=False
        image_tile = image.crop(tile)
        true_image_tile = image_tile.getbbox()
        if not true_image_tile:
            print "empty,skipping"
            continue
        tile_coords = [[true_image_tile[0], true_image_tile[1]],[true_image_tile[2], true_image_tile[1]],[true_image_tile[2], true_image_tile[3]], [true_image_tile[0], true_image_tile[3]]]
        tile_poly = shapely.geometry.MultiPoint([shapely.geometry.Point(x) for x in tile_coords]).convex_hull
        for label in labels:
            polygon_coords = [shapely.geometry.Point(x['x'],x['y']) for x in label['vertices']]
            polygon = shapely.geometry.MultiPoint(polygon_coords).convex_hull
            try:
                intersect = tile_poly.intersection(polygon)
            except:
                print "errored on intersect"
                continue
            #if type(intersect) == shapely.geometry.collection.GeometryCollection:
            if isinstance(intersect, shapely.geometry.collection.GeometryCollection) or isinstance(intersect, shapely.geometry.MultiPolygon):
                for entry in intersect:
                    print entry
                    print entry.area
            elif intersect.area > 30:
                if label["label_class"] == "land":
                    if (tile_poly.intersection(polygon).area/tile_poly.area)*100 > 80:
                        print "all_land"
                        #we don't want this much land, so we can skip this image
                        skip = True
                        break
                segmentations.append({
                        "category" : label["label_class"],
                        "segmentation" : [(x[0] - tile[0], x[1] - tile[1]) for x in numpy.array(intersect.exterior)],
                        "bbox" : [a - b for a, b in zip(intersect.bounds, [tile[0], tile[1], tile[0], tile[1]])],
                        "area" : intersect.area,
                        })
                """
                print "***************"
                print tile
                print "----------------"
                print tile_coords
                print tile_poly.bounds
                print intersect.bounds
                print [a - b for a, b in zip(intersect.bounds, [tile[0], tile[1], tile[0], tile[1]])]
                """
        if len(segmentations) and not skip:
            image_tile.save(os.path.join(output_location, os.path.basename(image_file.replace(".png", "_id{}.png".format(index+1)))))
            tile_json = {
                    "segmentations" : segmentations,
                    "image_id" : os.path.basename(image_file.replace(".png", "_id{}.png".format(index+1))),
                    "true_bounds" : tile
            }
            with open(os.path.join(output_location, os.path.basename(image_file.replace(".png", "_id{}.json".format(index+1)))), "w") as output:
                json.dump(tile_json, output)

            if spectral_data:
                sub_spectral_data = {}
                for band in spectral_data.keys():
                    if spectral_data[band].shape == (2000,2000):
                        sub_spectral_data[band] = spectral_data[band][tile[0]:tile[2], tile[1]:tile[3]]
                oname = os.path.join(output_location, os.path.basename(image_file.replace(".png", "_id{}.nc".format(index+1))))
                make_spectral_nc(sub_spectral_data, oname)
                
    return
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='create some sentinel subtiles')
    parser.add_argument('output', help="output location")
    parser.add_argument("image_files_list", help="text file of paths to files")
    parser.add_argument("--no_thread", action="store_true", help="text file of paths to files")
    args = parser.parse_args()
    print "json"
    if os.path.exists(os.path.join(args.output, "mega_json.json")):
        with open(os.path.join(args.output, "mega_json.json")) as jfile:
            mega_json = json.load(jfile)
    else:
        mega_json = {"annotations" : []}
    print "image"
    image_files_list = list(open(args.image_files_list))
    no = 1
    images = [x.replace("\n", "") for x in image_files_list]
    proc_args = [(x, args.output) for x in images]
    if not args.no_thread:
        pool = multiprocessing.Pool(4)
        pool.map(image_slicer, proc_args)
    else:
        for proc in proc_args:
            image_slicer(proc)
    """
    for image in image_files_list:
        print no
        no += 1
        image = image.replace("\n", "")
        t = multiprocessing.Process(target=image_slicer, args=(image, args.output))
        threads.append(t)
        image_slicer(image, args.output)
    """
    """
    print "json_2"
    for index, annotation in enumerate(mega_json["annotations"]):
        mega_json["annotations"][index]["id"] = index
    print "json_3"
    print mega_json
    with open(os.path.join(args.output, "mega_json.json"), "w") as output:
        json.dump(mega_json, output)
    print "finished!"
    """

