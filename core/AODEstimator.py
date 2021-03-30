# -*- coding: utf-8 -*-
# Class: AODEstimator
# Aeresol of Depth Estimation

import os
import math

import gdal
import numpy as np

import rasterio
import rasterio.plot
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from core.LandsatMetadataReader import LandsatMetadataReader
from core.CalibrateLandsatBand import CalibrateLandsatBand

class AODEstimator():

    def __init__(self, metadata_file, LSE_mode='auto-ndvi-raw', temp_dir=None):
        
        if LSE_mode not in self.LSE_modes:
            raise ValueError('Unsupported LSE mode. Supported: %s' % self.LSE_modes)

        self.metadata_file = metadata_file
        self.metadata = LandsatMetadataReader(self.metadata_file)

        self.dataset_basepath = os.path.dirname(self.metadata_file)

        self.band_10_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_10'])
        self.band_11_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_11'])
        
        print('Path to Band 10: ', self.band_10_path)
        print('Path to Band 11: ', self.band_11_path)

        self.temp_dir = temp_dir
        self.lse_mode = LSE_mode

    def transform_array_back_to_original_size(self, array, window_size, original_raster):
        original_geotransform = original_raster.GetGeoTransform()
        original_projection = original_raster.GetProjection()
        original_datatype = gdal.GDT_Float32

        xMin = original_geotransform[0]
        yMax = original_geotransform[3]
        xMax = xMin + original_geotransform[1] * original_raster.RasterXSize
        yMin = yMax + original_geotransform[5] * original_raster.RasterYSize

        new_geotransform = (original_geotransform[0],
                            original_geotransform[1] * window_size,
                            original_geotransform[2],
                            original_geotransform[3],
                            original_geotransform[4],
                            original_geotransform[5] * window_size)

        driver = gdal.GetDriverByName('MEM')
        new_ds = driver.Create('', array.shape[1], array.shape[0], 1, original_datatype)
        new_ds.GetRasterBand(1).WriteArray(array)
        new_ds.SetProjection(original_projection)
        new_ds.SetGeoTransform(new_geotransform)

        # to source grid
        new_ds_translated = gdal.Warp('', new_ds, format='MEM', xRes=original_geotransform[1],
                                      yRes=original_geotransform[5], outputBounds=[xMin, yMin, xMax, yMax],
                                      dstNodata=np.nan, resampleAlg='cubic')

        new_ds_translated_array = new_ds_translated.GetRasterBand(1).ReadAsArray()

        return new_ds_translated_array
        
    def __save_array_to_gtiff(self, array, domain_raster, gtiff_path):
        driver = gdal.GetDriverByName("GTiff")
        dataType = gdal.GDT_Float32
        dataset = driver.Create(gtiff_path, array.shape[1], array.shape[0], domain_raster.RasterCount, dataType)
        dataset.SetProjection(domain_raster.GetProjection())
        dataset.SetGeoTransform(domain_raster.GetGeoTransform())

        dataset.GetRasterBand(1).WriteArray(array)

        del dataset     