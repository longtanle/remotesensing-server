# -*- coding: utf-8 -*-
"""Handle login and downloading from the EarthData NASA portal."""
#!/usr/bin/env python
# get_modis A MODIS land product downloading tool
# Copyright (c) 2013-2016 J Gomez-Dans. All rights reserved.
#
# This file is part of get_modis.
#
# get_modis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# get_modis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with brdf_filter.  If not, see <http://www.gnu.org/licenses/>.
import optparse
import os
import os.path
import multiprocessing as mp

import subprocess, glob, shutil

import pyproj
import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

from osgeo import gdal, ogr, osr
from pyhdf.SD import SD, SDC

import rasterio
import rasterio as rst
from rasterio.merge import merge
from rasterio.plot import show

import time
import calendar
import logging
import sys
import fnmatch
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)


__author__ = "Long Tan Le"
__copyright__ = "Copyright 2020 Long Tan Le"
__version__ = "1.0.0"
__license__ = "GPLv3"
__email__ = "tanlong.ce@gmail.com"

LOG = logging.getLogger( __name__ )
OUT_HDLR = logging.StreamHandler( sys.stdout )
OUT_HDLR.setFormatter( logging.Formatter( '%(asctime)s %(message)s') )
OUT_HDLR.setLevel( logging.INFO )
LOG.addHandler( OUT_HDLR )
LOG.setLevel( logging.INFO )

HEADERS = { 'User-Agent' : 'earthdata Python %s' % __version__ }

CHUNKS = 65536

def return_url(url):
    the_day_today = time.asctime().split()[0]
    the_hour_now = int(time.asctime().split()[3].split(":")[0])
    if the_day_today == "Wed" and 14 <= the_hour_now <= 17:
        LOG.info("Sleeping for %d hours... Yawn!" % (18 - the_hour_now))
        time.sleep(60 * 60 * (18 - the_hour_now))

    req = urllib2.Request("%s" % (url), None, HEADERS)
    html = urllib2.urlopen(req).readlines()
    return html

def parse_modis_dates ( url, dates, product, out_dir, ruff=False ):
    """Parse returned MODIS dates.

    This function gets the dates listing for a given MODIS products, and
    extracts the dates for when data is available. Further, it crosses these
    dates with the required dates that the user has selected and returns the
    intersection. Additionally, if the `ruff` flag is set, we'll check for
    files that might already be present in the system and skip them. Note
    that if a file failed in downloading, it might still be around
    incomplete.

    Parameters
    ----------
    url: str
        A URL such as "http://e4ftl01.cr.usgs.gov/MOTA/MCD45A1.005/"
    dates: list
        A list of dates in the required format "YYYY.MM.DD"
    product: str
        The product name, MOD09GA.005
    out_dir: str
        The output dir
    ruff: bool
        Whether to check for present files
    Returns
    -------
    A (sorted) list with the dates that will be downloaded.
    """
    if ruff:
        product = product.split(".")[0]
        already_here = fnmatch.filter(os.listdir(out_dir),
                                      "%s*hdf" % product)
        already_here_dates = [x.split(".")[-5][1:]
                              for x in already_here]

    html = return_url(url)
    
    available_dates = []
    for line in html:

        if line.decode().find("href") >= 0 and \
                        line.decode().find("[DIR]") >= 0:
            # Points to a directory
            the_date = line.decode().split('href="')[1].split('"')[0].strip("/")
            if ruff:
                try:
                    modis_date = time.strftime("%Y%j",
                                               time.strptime(the_date,
                                                             "%Y.%m.%d"))
                except ValueError:
                    continue
                if modis_date in already_here_dates:
                    continue
                else:
                    available_dates.append(the_date)
            else:
                available_dates.append(the_date)

    dates = set(dates)
    available_dates = set(available_dates)
    suitable_dates = list(dates.intersection(available_dates))
    suitable_dates.sort()
    return suitable_dates

def get_modisfiles(username, password, product, year, tile, proxy,
                   doy_start=1, doy_end=-1,
                   base_url="https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD19A2/", out_dir=".",
                   ruff=False, get_xml=False, verbose=False):

    """Download MODIS products for a given tile, year & period of interest

    This function uses the `urllib2` module to download MODIS "granules" from
    the USGS website. The approach is based on downloading the index files for
    any date of interest, and parsing the HTML (rudimentary parsing!) to search
    for the relevant filename for the tile the user is interested in. This file
    is then downloaded in the directory specified by `out_dir`.

    The function also checks to see if the selected remote file exists locally.
    If it does, it checks that the remote and local file sizes are identical.
    If they are, file isn't downloaded, but if they are different, the remote
    file is downloaded.

    Parameters
    ----------
    username: str
        The EarthData username string
    password: str
        The EarthData username string
    platform: str
        One of three: MOLA, MOLT MOTA
    product: str
        The product name, such as MOD09GA.005 or MYD15A2.005. Note that you
        need to specify the collection number (005 in the examples)
    year: int
        The year of interest
    tile: str
        The tile (e.g., "h17v04")
    proxy: dict
        A proxy definition, such as {'http': 'http://127.0.0.1:8080', \
        'ftp': ''}, etc.
    doy_start: int
        The starting day of the year.
    doy_end: int
        The ending day of the year.
    base_url: str, url
        The URL to use. Shouldn't be changed, unless USGS change the server.
    out_dir: str
        The output directory. Will be create if it doesn't exist
    ruff: Boolean
        Check to see what files are already available and download them without
        testing for file size etc.
    verbose: Boolean
        Whether to sprout lots of text out or not.
    get_xml: Boolean
        Whether to get the XML metadata files or not. Someone uses them,
        apparently ;-)
    Returns
    -------
    Nothing
    """

    if proxy is not None:
        proxy = urllib2.ProxyHandler(proxy)
        opener = urllib2.build_opener(proxy)
        urllib2.install_opener(opener)

    if not os.path.exists(out_dir):
        if verbose:
            LOG.info("Creating outupt dir %s" % out_dir)
        os.makedirs(out_dir)
    if doy_end == -1:
        if calendar.isleap(year):
            doy_end = 367
        else:
            doy_end = 366

    dates = [time.strftime("%Y.%m.%d", time.strptime("%d/%d" % (i, year),
                                                     "%j/%Y")) for i in
             range(doy_start, doy_end)]
    # url = "%s/%s/%s/" % (base_url, platform, product)
    url = base_url
    print("URL here " + url)
    
    #dates = parse_modis_dates(url, dates, product, out_dir, ruff=ruff)

    print(dates)
    
    them_urls = []
    store_path = []
    
    r = requests.get(url, verify=False)
    for line in r.text.split("\n"):
        if line.find(tile) >= 0:
            if line.find(".hdf")  >= 0:
                fname = line.split("data-name=")[1].split(" ")[0].strip('"')
                print(fname)
                if not os.path.exists(os.path.join(out_dir, fname)):
                    them_urls.append("%s/%s" % (url, fname))
                else:
                    if verbose:
                        LOG.info("File %s already present. Skipping" % fname)
                store_path.append("%s/%s" % (out_dir,fname))
                break;
    print(them_urls)
    with requests.Session() as s:
        adapter = requests.adapters.HTTPAdapter(max_retries=3)
        s.mount('https://', adapter)
        s.mount('http://', adapter)
        s.auth = (username, password)
        for the_url in them_urls:
            r1 = s.request('get', the_url)
            try:
                r = s.get(r1.url, allow_redirects=True, stream=True, timeout=300)

                if not r.ok:
                    raise IOError("Can't start download... [%s]" % the_url)
                file_size = int(r.headers.get("Content-Length"))
                fname = the_url.split("/")[-1]
                LOG.info("Starting download on %s(%d bytes) ..." %
                        (os.path.join(out_dir, fname), file_size))
                with open(os.path.join(out_dir, fname), 'wb') as fp:
                    for chunk in r.iter_content(chunk_size=CHUNKS):
                        if chunk:
                            fp.write(chunk)
                    fp.flush()
                    os.fsync(fp)
                    if verbose:
                        LOG.info("\tDone!")
            except requests.exceptions.Timeout:
                LOG.info("Connection timeout after {} seconds.".format(timeout))
    if verbose:
        LOG.info("Completely finished downlading all there was")
        
    return store_path


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# CONVERSION STRATEGY --> MODIS. 
# Using GDAL CLI tools
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def extract_bands( fn ):
	print( fn )
	# translate to GTiff
	return subprocess.call([ 'gdal_translate', '-q', '-of', 'GTiff', '-sds', '-co', 'COMPRESS=LZW', \
								fn, os.path.basename(fn).replace('.','_').replace('_hdf','.tif') ])

def run_all_extract_bands( files, output_dir, ncpus=16 ):
	# change the directory
	os.chdir( output_dir )
	# parallel it
	pool = mp.Pool( ncpus )
	out = pool.map( extract_bands, files )
	pool.close()
	pool.join()
	return output_dir

def move_file( x ):
	fn, out_fn = x
	return shutil.move( fn, out_fn )

def move_desired_bands( in_dir, out_dir, keep_bands=['01','03'], ncpus=14 ):
	
	# list and filter the files based on the desired bands
	files = glob.glob( os.path.join( in_dir, '*.tif' ) )
	files = [fn for fn in files if fn.endswith(tuple('_{}.tif'.format(band) for band in keep_bands))]
	out_files = [ os.path.join( out_dir, os.path.basename(fn) ) for fn in files ]
	# move the files to a new dir
	pool = mp.Pool(ncpus)
	moved = pool.map(move_file, list(zip(files,out_files)))
	pool.close()
	pool.join()
	return out_dir

def move_undesired_bands( in_dir, out_dir, ncpus=16 ):
	files = glob.glob(os.path.join(in_dir,'*.tif'))
	out_files = [ os.path.join( out_dir, os.path.basename(fn) ) for fn in files ]
	pool = mp.Pool(ncpus)
	out = pool.map(move_file, list(zip(files, out_files)))
	pool.close()
	pool.join()
	return out_dir	

def make_band_lookup():
	'''hardwired name lookup table for M*D11A2 data...'''
	return {'01':'Optical_Depth_047','02':'Optical_Depth_055','03':'AOD_Uncertainty','04':'FineModeFraction',
			'05':'Column_WV','06':'AOD_QA','07':'AOD_MODEL','08':'Injection_Height',
			'09':'cosSZA','10':'cosVZA','11':'RelAZ','12':'Scattering_Angle','13':'Glint_Angle'}

def make_mosaic_args( files, out_dir ):
	'''make arguments to pass to mosaic'''

	# make a dataframe that we can group and use in arguments creation for the mosaicking
	colnames = ['product', 'date', 'tile', 'version', 'production', 'band']
	df = pd.DataFrame([os.path.basename(fn).split('.')[0].split('_') for fn in files], columns=colnames )
	df['fn'] = files

	def make_args(group, out_dir):
		files = sorted(group['fn'].tolist())
		elems = group[['product', 'date','version', 'production', 'band']].iloc[0].tolist()
		out_fn = os.path.join( out_dir, '_'.join(elems) + '.tif').replace('_006_','_HCM_006_')
		return [files] + [out_fn]

	return [make_args( group, out_dir ) for i,group in df.groupby([ 'product', 'date', 'band' ])]

def mosaic_tiles(files, out_fn):
	command = ['gdal_merge.py','-n','0','-a_nodata','0', '-o', out_fn,] + files
	_ = subprocess.call(command)
	return out_fn

def wrap_mosaic_tiles( x ):
	'''a wrapper for the mosaic f(x) for parallelism'''
	files, out_fn = x
	return mosaic_tiles( files, out_fn )

def run_mosaic_tiles( args, ncpus=5 ):
	pool = mp.Pool( ncpus )
	out = pool.map( wrap_mosaic_tiles, args )
	pool.close()
	pool.join()

def warp_to_3338( fn, out_fn ):
	return subprocess.call([ 'gdalwarp','-q','-overwrite',
				'-t_srs','EPSG:32648','-co','COMPRESS=LZW', fn, out_fn ])

def wrap_warp_to_3338(x):
	return warp_to_3338(*x)

def run_warp_to_3338( args, ncpus=5 ):
	pool = mp.Pool( ncpus )
	out = pool.map( wrap_warp_to_3338, args )
	pool.close()
	pool.join()
	return out

def rescale_values( fn ):

	print(fn)
	with rasterio.open( fn ) as rst:
		meta = rst.meta.copy()
		print(meta)
        
		meta.update( compress='lzw', dtype='float32', nodata=0 )
		arr = rst.read(1)

	print(arr.shape)
	arr_out = np.copy(arr).astype(np.double)

	# scale it:
	if fn.endswith(('_01.tif','_02.tif','_03.tif','_04.tif','_09.tif','_10.tif')):
		
	 	invalid = np.logical_or(arr_out < -100, arr_out > 5000)
	 	invalid = np.logical_or(invalid, arr_out == -28672)
	 	arr_out[invalid] = np.nan
	 	arr_out = arr_out * 0.001
	 	arr_out = np.ma.masked_array(arr_out, np.isnan(arr_out))

	 	print(arr_out)
	
	# elif fn.endswith(('_11.tif','_12.tif','_13.tif')):
	# 	invalid = np.logical_or(arr_out < -100, arr_out > 5000)
	# 	invalid = np.logical_or(invalid, arr_out == -28672)
	# 	arr_out[invalid] = np.nan
	# 	arr_out = arr_out * 0.01
	# 	arr_out = np.ma.masked_array(arr_out, np.isnan(arr_out))

	else:
		raise BaseException('wrong bands')

    

	# make the output filename and dump to disk
	out_fn = fn.replace( 'warped', 'rescaled' )
	with rasterio.open( out_fn, 'w', **meta ) as out:
	 	plt.imshow(arr_out, cmap='RdYlGn')
	 	plt.colorbar()
	 	plt.title('Rescaled')
	 	plt.xlabel('Column #')
	 	plt.ylabel('Row #')           
	 	png_fn = out_fn.replace( '.tif', '.png' )
	 	plt.savefig(png_fn)          
	 	out.write( arr_out.astype( np.float32 ), 1 )
	return out_fn


def run_rescale_values( files, ncpus ):
	pool = mp.Pool( ncpus )
	out = pool.map( rescale_values, files )
	pool.close()
	pool.join()
	return out

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   
def save_raster (raster_data, output_name, dataset, NaN_Value):
    '''save a np.array using GDAL to a GeoTiff file
    raster_data is a np.array
    output_name is the file to store
    dataset is a GDAL-accesible dataset to collect the projection and geotransform from
    NaN_Value defines what is placed
    '''
    # Open the reference dataset
    g = (dataset)
    # Get the Geotransform vector
    geo_transform = g.GetGeoTransform()
    x_size = g.RasterXSize # Raster xsize
    y_size = g.RasterYSize # Raster ysize
    srs = 'PROJCS["WGS 84 / Pseudo-Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]],PROJECTION["Mercator_1SP"],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"],AUTHORITY["EPSG","3857"]]'
    #srs = g.GetProjectionRef () # Projection
    #raster_data = g.ReadAsArray().astype(float)
    NaN_rast = NaN_Value
    # raster_data[raster_data == NaN_rast] = 'NaN'
    raster_data[raster_data == NaN_rast] = np.nan
    # Need a driver object. By default, we use GeoTIFF
    dst_options = ['COMPRESS=LZW']
    driver = gdal.GetDriverByName("GTiff")
    dataset_out = driver.Create(output_name, x_size, y_size, 1, gdal.GDT_Float32, dst_options)
    dataset_out.SetGeoTransform(geo_transform)
    dataset_out.SetProjection(srs)
    dataset_out.GetRasterBand(1).WriteArray(raster_data.astype(np.float32))
    dataset_out.FlushCache()
    
    return

def save_raster_comp (raster_data, bands, output_name, dataset, NaN_Value):
    '''save a np.array using GDAL to a GeoTiff file
    raster_data is a np.array
    output_name is the file to store
    dataset is a GDAL-accesible dataset to collect the projection and geotransform from
    NaN_Value defines what is placed
    '''
    # Open the reference dataset
    g = (dataset)
    # Get the Geotransform vector
    geo_transform = g.GetGeoTransform()
    x_size = g.RasterXSize # Raster xsize
    y_size = g.RasterYSize # Raster ysize
    
    print(x_size)
    print(y_size)
    
    srs = 'PROJCS["WGS_1984_UTM_Zone_48N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
    #srs = g.GetProjectionRef () # Projection
    #raster_data = g.ReadAsArray().astype(float)
    NaN_rast = NaN_Value
    # raster_data[raster_data == NaN_rast] = 'NaN'
    # Need a driver object. By default, we use GeoTIFF
    dst_options = ['COMPRESS=LZW']
    driver = gdal.GetDriverByName("GTiff")
    dataset_out = driver.Create(output_name, x_size, y_size, bands, gdal.GDT_Float32, dst_options)
    dataset_out.SetGeoTransform(geo_transform)
    dataset_out.SetProjection(srs)
    for b in np.arange(bands):
        rst = raster_data[:,:,b]
        rst[rst == NaN_rast] = np.nan
        dataset_out.GetRasterBand(b+1).WriteArray(rst.astype(np.float32))
    dataset_out.FlushCache()
    
    return

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# read and process set of terra and aqua tiles
# store to GeoTiff
# platform, product, year, tile, proxy,
#                   doy_start=1, doy_end=-1,
#                   base_url="http://e4ftl01.cr.usgs.gov", out_dir=".",
#                   ruff=False, get_xml=False, verbose=False
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def tile_process(store_path_terr,store_path_aqua,ti):
    '''Processing of one tile over one month
    fi_terra and fi_aqua are the respective file arrays
    ti is the tile index
    
    This function calls to download the hdf files. 
    It then extracts the required LST data and QC indicees.
    A np.array is filled with the requested data.
    The mean values for day and night are stored to two GeoTiff files.'''
    
    #store_path_terr = download_fi(fi_terra[ti])
    #store_path_aqua = download_fi(fi_aqua[ti])
    
    # create empty arrays
    dt_arr = np.zeros((1200,1200,len(store_path_terr)))
    nt_arr = np.zeros((1200,1200,len(store_path_terr)))
    
    da_arr = np.zeros((1200,1200,len(store_path_terr)))
    na_arr = np.zeros((1200,1200,len(store_path_terr)))


    # get data from terrestial LST
    for i in np.arange(len(store_path_terr)):
        
        try:
            # load layers from hdf4 file
            gdal_dataset = gdal.Open(store_path_terr[i])
            subsets = gdal_dataset.GetSubDatasets()
            dataset_LSTd = gdal.Open(subsets[0][0])
            qc_LSTd = gdal.Open(subsets[1][0])
            dataset_LSTn = gdal.Open(subsets[4][0])
            qc_LSTn = gdal.Open(subsets[5][0])
        
            # populate array with data wgridmetahere the criterion (QC = 0) is met
            dummy = dataset_LSTd.ReadAsArray().astype(float)
            
            dummy[qc_LSTd.ReadAsArray()>0] = np.nan
            dt_arr[:,:,i] = dummy
            dummy = dataset_LSTn.ReadAsArray().astype(float)
            dummy[qc_LSTn.ReadAsArray()>0] = np.nan
            nt_arr[:,:,i] = dummy
            
        except:
            print('Error while reading '+ store_path_terr[i] + '.')
            print('Tile data not considered.')
            dt_arr[:,:,i] = np.nan
            nt_arr[:,:,i] = np.nan

    # get data from aquatic LST
    for i in np.arange(len(store_path_aqua)):
        
        try:
            # load layers from hdf4 file
            gdal_dataset = gdal.Open(store_path_aqua[i])
            subsets = gdal_dataset.GetSubDatasets()
            dataset_LSTd = gdal.Open(subsets[0][0])
            qc_LSTd = gdal.Open(subsets[1][0])
            dataset_LSTn = gdal.Open(subsets[4][0])
            qc_LSTn = gdal.Open(subsets[5][0])
            
            # populate array with data where the criterion (QC = 0) is met
            dummy = dataset_LSTd.ReadAsArray().astype(float)
            dummy[qc_LSTd.ReadAsArray()>0] = np.nan
            da_arr[:,:,i] = dummy
            dummy = dataset_LSTn.ReadAsArray().astype(float)
            dummy[qc_LSTn.ReadAsArray()>0] = np.nan
            na_arr[:,:,i] = dummy
        except:
            print('Error while reading ' + store_path_aqua[i] + '.')
            print('Tile data not considered.')
            da_arr[:,:,i] = np.nan
            na_arr[:,:,i] = np.nan

    # calculate mean of LST records in time window
    d_arr_out = np.zeros((1200,1200,2))
    n_arr_out = np.zeros((1200,1200,2))
    # Compute the arithmetic mean along the specified axis, ignoring NaNs.
    d_arr_out[:,:,0] = np.nanmean(dt_arr,axis=2)
    d_arr_out[:,:,1] = np.nanmean(da_arr,axis=2)
    n_arr_out[:,:,0] = np.nanmean(nt_arr,axis=2)
    n_arr_out[:,:,1] = np.nanmean(na_arr,axis=2)

    #save_raster (d_arr_out, 2, 'res/'+str(ti).zfill(4)+'_d.tif', dataset_LSTd, 0)
    #save_raster (n_arr_out, 2, 'res/'+str(ti).zfill(4)+'_n.tif', dataset_LSTn, 0)
    
    #save_raster_comp (d_arr_out,2, 'res/' + str(ti).zfill(4) +'_d.tif', dataset_LSTd, 0)
    #save_raster_comp (n_arr_out,2, 'res/' + str(ti).zfill(4) + '_n.tif', dataset_LSTn, 0)
    save_raster_comp (dt_arr,2, 'res/' + str(ti).zfill(4) +'_d.tif', dataset_LSTd, 0)
    save_raster_comp (nt_arr,2, 'res/' + str(ti).zfill(4) + '_n.tif', dataset_LSTn, 0)
    
    return ['res/'+ str(ti).zfill(4) +'_d.tif', 'res/'+ str(ti).zfill(4) +'_n.tif']

def mosaic_Gtif(tis,out_fi='mos.tif'):
    '''Mosaic GeoTiff files
    tis is a list of GeoTiff files to be stitched together
    out_fi is the name of the output GeoTiff file
    '''
    
    mos_fi = tis#['test1.tif','test2.tif']
    fi_to_mosaic = []
    for fi in mos_fi:
        sx = rst.open(fi)
        fi_to_mosaic.append(sx)
    
    mosaic, out_trans = merge(fi_to_mosaic)
    show(mosaic)

    # write mosaic to file
    out_meta = sx.meta.copy()
    out_meta.update({"driver": "GTiff","height": mosaic.shape[1],"width": mosaic.shape[2],"transform": out_trans,"compress": "LZW","crs": "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"})
    #'crs': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    #"crs": "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs "
    
    with rst.open(out_fi, "w", **out_meta) as dest:
        dest.write(mosaic)
        
    return


def wrap_tile_process(store_path_terr, sgridmetatore_path_aqua, ti):
    #if os.path.isfile('./res/'+str(ti).zfill(4)+'_d.tif')==False:
    try:
        tile_process(store_path_terr,store_path_aqua,ti)
    except:
        print(ti+' could not be processed')
    return

def processHDF(filename):
    # Read HDF File 
    hdf = SD(filename, SDC.READ)

    datasets = hdf.datasets()

    for i,v in enumerate(datasets):
        print('{0}. {1}'. format(i+1,v))
    
    # Select Optical Depth .55 micron wavelength
    optical_depth_055 = hdf.select('Optical_Depth_055')
    optical_depth = optical_depth_055[0,:,:].astype(np.double)
    shape = optical_depth_055[:,:,:].shape

    # Obtain fillValue and Scale factor for optical depth
    OD_fillValue  = optical_depth_055.getfillvalue()
    OD_scaleFactor = optical_depth_055.attributes(full=1)["scale_factor"][0]
    
    attrs = optical_depth_055.attributes(full=1)
    lna=attrs["long_name"]
    long_name = lna[0]
    vra=attrs["valid_range"]
    valid_range = vra[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]        
    ua=attrs["unit"]
    units = ua[0]
    aoa=attrs["add_offset"]
    add_offset = aoa[0]

    # Flatten and merge data based on number of orbits
    # optical_depth = np.array([], dtype=np.int16).reshape(0,1)
    # for i in range(shape[0]):
    #    optical_depth = np.append(optical_depth,optical_depth_055[i,:,:].ravel())

    # Scale data based on the scale factor.

    # Apply the attributes to the data.
    invalid = np.logical_or(optical_depth < valid_range[0], optical_depth > valid_range[1])
    invalid = np.logical_or(invalid, optical_depth == _FillValue)
    optical_depth[invalid] = np.nan
    optical_depth = (optical_depth - add_offset) * scale_factor
    optical_depth = np.ma.masked_array(optical_depth, np.isnan(optical_depth))

    print(OD_fillValue)
    print(OD_scaleFactor)
    print(optical_depth)

    # Select water vopur readings
    Column_WV_n = hdf.select('Column_WV')

    # Obtain fillValue and Scale factor for optical depth
    WV_fillValue  = Column_WV_n.getfillvalue()
    WV_scaleFactor = Column_WV_n.attributes(full=1)["scale_factor"][0]

    # Flatten and merge data based on number of orbits
    Column_WV = np.array([], dtype=np.int16).reshape(0,1)
    for i in range(shape[0]):
        Column_WV = np.append(Column_WV,Column_WV_n[i,:,:].ravel())
    
    # Scale data based on the scale factor.
    Column_WV = Column_WV*WV_scaleFactor

    # Select water vopur readings
    AOD_Uncertainty_n = hdf.select('AOD_Uncertainty')

    # Obtain fillValue and Scale factor for optical depth
    UN_fillValue  = AOD_Uncertainty_n.getfillvalue()
    UN_scaleFactor = AOD_Uncertainty_n.attributes(full=1)["scale_factor"][0]

    # Flatten and merge data based on number of orbits
    AOD_Uncertainty = np.array([], dtype=np.int16).reshape(0,1)
    for i in range(shape[0]):
        AOD_Uncertainty = np.append(AOD_Uncertainty,AOD_Uncertainty_n[i,:,:].ravel())

    # Scale data based on the scale factor.
    AOD_Uncertainty = AOD_Uncertainty*UN_scaleFactor
    
    return "TestModisAOD"
