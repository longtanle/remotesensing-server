# -*- coding: utf-8 -*-
# Class: LSTEstimator
# Land Surface Temperature Estimation

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


class LSTEstimator():
    
    # NDVI_max & NDVI_min
    NDVIs = 0.2
    NDVIv = 0.5
    
    # The emissivity of vegetation and soil
    e_v = 0.985
    e_s = 0.960
    
    # Wavelengths
    w_10 = 1.090e-5
    w_11 = 1.205e-5
    # Hang so Stefan - Boltzmann
    o = 1.380649e-23
    # Hang so Plank
    h = 6.62607004e-34
    # Toc do anh sang
    c = 2.998e8
    
    p = (h * c) / o
    
    # For calculating cloud & water vapor

    # Tinh LSE tu NDVI
    LSE_modes = ['auto-ndvi-raw']
    
    DEBUG = True

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

        self.scene_id = self.metadata.metadata['LANDSAT_SCENE_ID']
        self.temp_dir = temp_dir
        self.lse_mode = LSE_mode


    # Tinh LSE tu NDVI va Pv
    # 
    #
    def get_lse_from_ndvi(self, srem=False):
        self.band_4_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_4'])
        self.band_5_path = os.path.join(self.dataset_basepath, self.metadata.metadata['FILE_NAME_BAND_5'])

        print("Calibrating Band 4...")
        band4_calibrator = CalibrateLandsatBand(self.band_4_path, self.metadata_file)
        band4_reflectance = band4_calibrator.get_reflectance_as_array_2()
    #    band4_reflectance[band4_reflectance > 1] = 1
        
        print("Calibrating Band 5...")
        band5_calibrator = CalibrateLandsatBand(self.band_5_path, self.metadata_file)      
        band5_reflectance = band5_calibrator.get_reflectance_as_array_2()
    #    band5_reflectance[band5_reflectance > 1] = 1
        
        if self.DEBUG == True:
            self.__save_array_to_gtiff(band4_reflectance, gdal.Open(self.band_10_path),
                                       os.path.join(self.temp_dir, self.scene_id + '_Reflectance_B4.tif'))
            self.__save_array_to_gtiff(band5_reflectance, gdal.Open(self.band_10_path),
                                       os.path.join(self.temp_dir, self.scene_id + '_Reflectance_B5.tif'))
        
        # NDIV = (NIR - RED)/(NIR + RED)
        ndvi = (band5_reflectance - band4_reflectance) / (band5_reflectance + band4_reflectance)
        
        print('Show NDVI array:')
        print('Size: ' + str(len(ndvi)) + 'x' + str(len(ndvi[0])) )

        
        # from Land Surface Temperature Retrieval from Landsat 5, 7, and 8 over Rural Areas: 
        # Assessment of Different Retrieval Algorithms and Emissivity Models and Toolbox Implementation
        
        # Tinh Vegetation Proportion (Pv)
        # Pv = [(NDVI - NDVImin)/(NDVImax - NDVImin)]^2
        Pv = ((ndvi - self.NDVIs) / (self.NDVIv - self.NDVIs)) ** 2
        
        print('Show Pv array:')
        print('Size: ' + str(len(Pv)) + 'x' + str(len(Pv[0])) )

        
        # Tinh LSE 
        # dE_b10 = (1 - 0.9668) * (1 - Pv) * 0.55 * 0.9863
        # dE_b11 = (1 - 0.9747) * (1 - Pv) * 0.55 * 0.9896
        
        lse = self.e_s * (1 - Pv) + self.e_v * Pv
        
        print('Show LSE array:')
        print('Size: ' + str(len(lse)) + 'x' + str(len(lse[0])) )

        
        if self.DEBUG == True:
            self.__save_array_to_gtiff(ndvi, gdal.Open(self.band_10_path), 
                                       os.path.join(self.temp_dir, self.scene_id + '_NDVI.tif'))

        return lse



    def get_b10_b11_brightness_temp_arrays(self):
        print("Calibrating Band 10...")
        b10_calibrator = CalibrateLandsatBand(self.band_10_path, self.metadata_file)
        b10_bt = b10_calibrator.get_brightness_temperature_as_array()
        
        print("Calibrating Band 11...")
        b11_calibrator = CalibrateLandsatBand(self.band_11_path, self.metadata_file) 
        b11_bt = b11_calibrator.get_brightness_temperature_as_array()


        if self.DEBUG  == True:
            self.__save_array_to_gtiff(b10_bt, gdal.Open(self.band_10_path),
                                       os.path.join(self.temp_dir, self.scene_id + '_BT_B10.tif'))
            self.__save_array_to_gtiff(b11_bt, gdal.Open(self.band_10_path),
                                       os.path.join(self.temp_dir, self.scene_id + '_BT_B11.tif'))

        return b10_bt, b11_bt

    def get_lst_array(self):
        
        print('Start calculating LST..')
        
        print('Calculating Brightness Temperatures...')
        b10_bt, b11_bt = self.get_b10_b11_brightness_temp_arrays()
        
        # LSE
        
        print('Calculating LSE...')
        if self.lse_mode == 'auto-ndvi-raw':
            lse = self.get_lse_from_ndvi(srem=False)
        
        
        print('Estimating LST...')
        
        print(self.p)
        
        # CALCULATE LAND SURFACE TEMPERATURE
        lst_b10 = b10_bt / (1 + ((self.w_10 * b10_bt) / self.p ) * np.log(lse))
        lst_b11 = b11_bt / (1 + ((self.w_11 * b11_bt) / self.p ) * (np.log(lse)))    
        
        print('Show LST Band 10:')
        print('LST MIN VLUES : ' + str(np.amin(lst_b10)))
        print('LST MAX VLUES : ' + str(np.amax(lst_b10)))
        
        # plt.imshow(lst_b10, cmap='RdYlGn')
        # plt.colorbar()
        # plt.title('LST B10')
        # plt.xlabel('Column #')
        # plt.ylabel('Row #')
        # plt.savefig(os.path.join(self.temp_dir, 'lst_b10.png'))
        
        # plt.imshow(lst_b11, cmap='RdYlGn')
        # plt.colorbar()
        # plt.title('LST B11')
        # plt.xlabel('Column #')
        # plt.ylabel('Row #')
        # plt.savefig(os.path.join(self.temp_dir, 'lst_b11.png'))
        
        if self.DEBUG  == True:
            self.__save_array_to_gtiff(lse, gdal.Open(self.band_10_path),
                                       os.path.join(self.temp_dir, self.scene_id + '_LSE.tif'))
            self.__save_array_to_gtiff(lst_b10, gdal.Open(self.band_10_path),
                                           os.path.join(self.temp_dir, self.scene_id + '_LST_B10.tif'))  
            self.__save_array_to_gtiff(lst_b11, gdal.Open(self.band_11_path),
                                           os.path.join(self.temp_dir, self.scene_id + '_LST_B11.tif'))
        
        return lst_b10, lst_b11

    def get_lst_as_gtiff(self, output_path):
        lst = self.get_lst_array()
        self.__save_array_to_gtiff(lst, gdal.Open(self.band_10_path), os.path.join(self.temp_dir, output_path))

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

        # to source grid dstNodata=np.nan
        new_ds_translated = gdal.Warp('', new_ds, format='MEM', xRes=original_geotransform[1],
                                      yRes=original_geotransform[5], outputBounds=[xMin, yMin, xMax, yMax],
                                      dstNodata=-99, resampleAlg='cubic')

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
