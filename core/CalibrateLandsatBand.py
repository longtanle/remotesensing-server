# -*- coding: utf-8 -*-
# Class: CalibrateLandsatBand
# Hieu chinh cac Band trong Dataset
import os
import math

# Su dung Gdal de xu ly anh Landsat
import gdal
import numpy as np

from core.LandsatMetadataReader import LandsatMetadataReader

class CalibrateLandsatBand():
    # Khoi tao, doc file Metadata va cac Band can xu ly
    
    DEBUG = False
    
    def __init__(self, band_file, metadata_file):
        self.metadata_reader = LandsatMetadataReader(metadata_file)

        self.metadata = self.metadata_reader.metadata

        self.band_metadata = self.metadata_reader.get_band_metadata_by_file_name(band_file)

        if not self.band_metadata:
            raise KeyError('Invalid band')


        self.band_dataset = gdal.Open(band_file)
        self.band_array = self.band_dataset.GetRasterBand(1).ReadAsArray()
        
        print('Size: ' + str(len(self.band_array)) + 'x' + str(len(self.band_array[0])) )
        print('QCAL MIN VLUES : ' + str(np.amin(self.band_array)))
        print('QCAL MAX VLUES : ' + str(np.amax(self.band_array)))

    # HIEU CHINH BUC XA OLI BAND (Radiometric Calibratation)
    # Phuong phap: Gain and Bias
    # L_i = M_i * Q_CAL + A_i
    def get_radiance_as_array_2(self):
        if self.DEBUG == True:        
            print('Reflectance Mult: ', self.band_metadata['reflectance_mult'])
            print('Reflectance Add: ', self.band_metadata['reflectance_add'])
        #    print('Quantize Cal: ', self.band_metadata['quantize_cal_maximum'])
        #    print('Band Array: ', self.band_array)
        
        radiance = (self.band_metadata['reflectance_mult']) * (self.band_array) + self.band_metadata['reflectance_add']
        radiance[self.band_array==0] = np.nan
        return radiance
    
    # Radiance to ToA Reflectence
    # Phuong phap 1
    # p_i = pi * L_i * d^2 / ESUN_i * cos(solar_zenith_angle)
    def get_reflectance_as_array(self, not_native_radiance_array=False):
        if self.band_metadata['type'] != 'reflectance':
            raise TypeError('Given band is thermal')
        if type(not_native_radiance_array)==bool:
            radiance = self.get_radiance_as_array()
        else:
            radiance = not_native_radiance_array
        d = float(self.metadata['EARTH_SUN_DISTANCE'])
        O = np.deg2rad(float(self.metadata['SUN_ELEVATION']))
        E = self.band_metadata['solar_irradiance']

        reflectance = (np.pi*radiance*d*d)/(E*np.sin(O))
        return reflectance

    # Radiance to ToA Reflectence
    # Phuong phap 2
    # p_i = p'_i/sin(SE)
    def get_reflectance_as_array_2(self, not_native_radiance_array=False):
        if self.band_metadata['type'] != 'reflectance':
            raise TypeError('Given band is thermal')
        if type(not_native_radiance_array)==bool:
            radiance = self.get_radiance_as_array_2()
        else:
            radiance = not_native_radiance_array
        O = np.deg2rad(float(self.metadata['SUN_ELEVATION']))
        
        print('SUN ELEVATION:' + str(O))

        reflectance = radiance/(np.sin(O))
        return reflectance

    # HIEU CHINH BUC XA TIRS BAND (ATMOSPHERIC CORRECTION)
    # Phuong phap: Gain and Bias
    # L_i = M_i * Q_CAL + A_i
    def get_radiance_as_array(self):
        
        if self.DEBUG == True:
            print('Radiance Maximum: ', self.band_metadata['radiance_maximum'])
            print('Radiance Maximum: ', self.band_metadata['radiance_minimum'])
        #    print('Quantize Cal Maximum: ', self.band_metadata['quantize_cal_maximum'])
        #    print('Quantize Cal Minimum: ', self.band_metadata['quantize_cal_minimum'])
        #    print('Band Array: ', self.band_array)
            
    #   radiance = ((self.band_metadata['radiance_maximum']-self.band_metadata['radiance_minimum']) / (self.band_metadata['quantize_cal_maximum']-self.band_metadata['quantize_cal_minimum'])) * (self.band_array - self.band_metadata['quantize_cal_minimum']) + self.band_metadata['radiance_minimum']
        radiance = (self.band_metadata['radiance_mult']) * (self.band_array) + self.band_metadata['radiance_add']
        radiance[self.band_array==0] = np.nan
        return radiance
    
    # TINH NHIET DO SANG (BRIGHTNESS TEMPERATURE)
    def get_brightness_temperature_as_array(self):
        if self.band_metadata['type'] != 'thermal':
            raise TypeError('Given band is reflectance')

        radiance = self.get_radiance_as_array()

        K1 = self.band_metadata['k1_constant']

        K2 = self.band_metadata['k2_constant']
        
        print('K1 constant:' + str(K1))    
        print('K2 constant:' + str(K2))
            
        brightness_temperature = (K2 / (np.log(K1/(radiance+1))))
        
        brightness_temperature_C = brightness_temperature - 273
        

        print('Show Brightness Temperature Array:')
        print('Size: ' + str(len(brightness_temperature_C)) + 'x' + str(len(brightness_temperature_C[0])) )
        print('NDVI MIN VLUES : ' + str(np.amin(brightness_temperature_C)))
        print('NDVI MAX VLUES : ' + str(np.amax(brightness_temperature_C)))
        
       #     print(brightness_temperature)
       #     print(brightness_temperature_C)
        
        return brightness_temperature_C

    # Luu anh sau xu ly duoi dinh dang .tiff
    def save_array_as_gtiff(self, array, new_file_path):
        driver = gdal.GetDriverByName("GTiff")
        dataType = gdal.GDT_Float32
        dataset = driver.Create(new_file_path, self.band_dataset.RasterXSize, self.band_dataset.RasterYSize, self.band_dataset.RasterCount, dataType)
        dataset.SetProjection(self.band_dataset.GetProjection())
        dataset.SetGeoTransform(self.band_dataset.GetGeoTransform())
        dataset.GetRasterBand(1).WriteArray(array)
        del dataset
