# -*- coding: utf-8 -*-
# Class: LandsatMetadataReader
# XU LY VA DOC DU LIEU TU FILE METADATA
#

import os
#import math

#import gdal
#import numpy as np

class LandsatMetadataReader():
    
    # KHOI TAO, LAY DUONG DAN DEN FILE METADATA
    def __init__(self, metadata_file_path):
        self.metadata_file_path = metadata_file_path
        self.metadata_file = open(self.metadata_file_path,'r')
        self.metadata = {}
        self.bands = {}
        
        # Khoi tao gia tri Extra-terrestrial Solar Irradiation (ESUN)
        solar_irradiances = {
                             'LANDSAT_8':
                                 {'1': 1895.33,
                                  '2': 2004.57,
                                  '3': 1820.75,
                                  '4': 1549.49,
                                  '5': 951.76,
                                  '6': 366.97,
                                  '7': 247.55,
                                  '8': 85.46,
                                  '9': 1723.88}
                             }
        # Khoi tao gia tri buoc song cua cac Band W = (Wmin + Wmax)/2
        wavelengths = {'LANDSAT_8':
                                 {'1': 0.44,
                                  '2': 0.48,
                                  '3': 0.56,
                                  '4': 0.655,
                                  '5': 0.865,
                                  '6': 1.61,
                                  '7': 2.20,
                                  '8': 0.59,
                                  '9': 1.37,
                                  '10': 10.9,
                                  '11': 12.05}
                        }

        # Khoi tao gia tri cua sensor
        center_sensor_zenith = {'LANDSAT_8':0.001}
        
        # Doc va luu tru Metadata
        for line in self.metadata_file.readlines():
            if (line.find('GROUP') >= 0) or (line.find('=') == -1):
                continue
            else:
                line_normalized = line.replace(' ','')
                items = line_normalized.split('=')
                self.metadata[items[0]] = items[1].replace('\n','').replace('\"','')
                
        # Kiem tra du lieu Landsat 8
        if not 'SPACECRAFT_ID' in self.metadata:
            raise KeyError('Invalid metadata file')
    
        all_bands = []
        reflectance_bands = []
        thermal_bands = []
        
        # 11 band: 9 band OLI, 2 band TIRS
        if self.metadata['SPACECRAFT_ID'] == 'LANDSAT_8':
            all_bands = [1,2,3,4,5,6,7,8,9,10,11]
            reflectance_bands = [1,2,3,4,5,6,7,8,9]
            thermal_bands = [10,11]

        if not all_bands:
            raise KeyError('Invalid metadata file')
        # Luu cac tham so cua tung Band
        for band in all_bands:
            self.bands[str(band)] = {}

            self.bands[str(band)]['file_name'] = self.metadata['FILE_NAME_BAND_%s' % band]
            self.bands[str(band)]['number'] = band
            self.bands[str(band)]['radiance_maximum'] = float(self.metadata['RADIANCE_MAXIMUM_BAND_%s' % band])
            self.bands[str(band)]['radiance_minimum'] = float(self.metadata['RADIANCE_MINIMUM_BAND_%s' % band])
            self.bands[str(band)]['quantize_cal_maximum'] = float(self.metadata['QUANTIZE_CAL_MAX_BAND_%s' % band])
            self.bands[str(band)]['quantize_cal_minimum'] = float(self.metadata['QUANTIZE_CAL_MIN_BAND_%s' % band])
            self.bands[str(band)]['radiance_mult'] = float(self.metadata['RADIANCE_MULT_BAND_%s' % band])
            self.bands[str(band)]['radiance_add'] = float(self.metadata['RADIANCE_ADD_BAND_%s' % band])

            self.bands[str(band)]['wavelength'] = wavelengths[self.metadata['SPACECRAFT_ID']][str(band)]

            self.bands[str(band)]['center_sensor_zenith'] = center_sensor_zenith[self.metadata['SPACECRAFT_ID']]


            if band in reflectance_bands:
                self.bands[str(band)]['saturation'] = self.metadata['SATURATION_BAND_%s' % band]
                self.bands[str(band)]['reflectance_maximum'] = float(self.metadata['REFLECTANCE_MAXIMUM_BAND_%s' % band])
                self.bands[str(band)]['reflectance_minimum'] = float(self.metadata['REFLECTANCE_MINIMUM_BAND_%s' % band])
                self.bands[str(band)]['reflectance_mult'] = float(self.metadata['REFLECTANCE_MULT_BAND_%s' % band])
                self.bands[str(band)]['reflectance_add'] = float(self.metadata['REFLECTANCE_ADD_BAND_%s' % band])

                self.bands[str(band)]['solar_irradiance'] = solar_irradiances[self.metadata['SPACECRAFT_ID']][str(band)]

                self.bands[str(band)]['type'] = 'reflectance'

            if band in thermal_bands:
                self.bands[str(band)]['k1_constant'] = float(self.metadata['K1_CONSTANT_BAND_%s' % band])
                self.bands[str(band)]['k2_constant'] = float(self.metadata['K2_CONSTANT_BAND_%s' % band])

                self.bands[str(band)]['type'] = 'thermal'
    
    # Doc file cac Band trong Dataset
    def get_band_metadata_by_file_name(self, file_name):
        for band in self.bands.keys():
            if os.path.basename(file_name) == self.bands[band]['file_name']:
                return self.bands[band]
