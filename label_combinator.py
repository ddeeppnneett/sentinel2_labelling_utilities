#!/usr/bin/env python
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

class labels_data(object):
    def __init__(self,label_json_path):
        self.underlying_json_txt = open(label_json_path).read()
        self.underlying_json = json.loads(self.underlying_json_txt)
        self.labels = self.underlying_json["labels"]
    
    def combine_identical_labels(self):
        return True
    
    def remove_impossible_labels(self):
        """
        Will strip down anything that can't exist on top of each other (e.g. water can't exist beside cloud, algae can't exist on land etc)

        """
        return True

    def 