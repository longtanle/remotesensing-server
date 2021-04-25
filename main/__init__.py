import datetime
import logging as log

import coloredlogs
from flask import Flask, render_template
from flask_bcrypt import Bcrypt
from flask_cors import CORS
from flask_jwt_extended import (
    JWTManager,
    create_access_token,
    get_jwt_identity,
    jwt_required,
)
from flask import send_file

import optparse
import os
import os.path
import multiprocessing as mp

import subprocess, glob, shutil

import subprocess, glob, shutil

import json
from core.ee_apis import API
import tempfile
from core.ee_earthexplorer import EarthExplorer, EarthExplorerError
import core.earthdata as earthdata
import core.earthdata_aod as earthdata_aod
from core.LSTEstimator import LSTEstimator

from osgeo import gdal, ogr, osr

import zipfile
import tarfile



# local imports
from config import app_config

from .db import MongoDB

def copy_file( x ):
    fn, out_fn = x

    return shutil.copyfile( fn, out_fn )

def copy_desired_tiff( in_dir, out_dir, keep_bands=[ '_LST_B10', '_LST_B11'], ncpus=14 ):
    # list and filter the files based on the desired bands
    files = glob.glob( os.path.join( in_dir, '*.tif' ) )
    files = [fn for fn in files if fn.endswith(tuple('{}.tif'.format(band) for band in keep_bands))]
    print(files)
    out_files = [ os.path.join( out_dir, os.path.basename(fn) ) for fn in files ]
    # move the files to a new dir
    pool = mp.Pool(ncpus)
    copied = pool.map(copy_file, list(zip(files,out_files)))
    pool.close()
    pool.join()
    
    return out_dir


def create_app(config_name):
    config_name = "dev" if not config_name else config_name
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(app_config[config_name])
    app.config.from_pyfile("config.py")
    app.config["log"] = log

    cors = CORS(app)
    app.config["CORS_HEADERS"] = "Content-Type"

    app.config["JWT_SECRET_KEY"] = "9MZbGqQHaC47SSKyKaTK"
    app.config["JWT_BLACKLIST_ENABLED"] = True
    app.config["JWT_BLACKLIST_TOKEN_CHECKS"] = ["access", "refresh"]
    app.config["JWT_ACCESS_TOKEN_EXPIRES"] = datetime.timedelta(days=10)
    app.config["jwt"] = JWTManager(app)
    app.config["flask_bcrypt"] = Bcrypt(app)
    jwt = app.config["jwt"]


    # Swagger UI config
    app.config.SWAGGER_UI_JSONEDITOR = True
    app.config.SWAGGER_UI_DOC_EXPANSION = "none"  # none, list, full

    with app.app_context():
        db = MongoDB()

    @jwt.token_in_blacklist_loader
    def check_if_token_in_blacklist(decrypted_token):
        blacklist = set()
        jti = decrypted_token["jti"]
        return jti in blacklist

    @app.route("/")
    def hello_world():
        # render home template
        return render_template('index.html')

    @app.route('/ee/search/<cloudcover>/<startDate>/<endDate>')
    def search(cloudcover, startDate,endDate):
        # Initialize a new API instance and get an access key

        username = app.config["USGS_USERNAME"]
        password = app.config["USGS_PASSWORD"]
        api = API(username, password)

        # Search for Landsat TM scenes
        scenes = api.search(
            dataset= 'landsat_8_c1',
            latitude= 10.7432,
            longitude= 106.875,
            start_date= startDate,
            end_date= endDate,
            max_cloud_cover= int(cloudcover)
        )

        #print (scenes)
        #print(f"{len(scenes)} scenes found.")
        response = []
        # Process the result
        for scene in scenes:
            item = {
                "acquisition_date": scene['acquisition_date'].strftime('%Y-%m-%d'),
                "scene_id": scene['entity_id'],
                "display_id": scene['display_id'],
                "wrs_path" : scene['wrs_path'],
                "wrs_row" : scene['wrs_row'] 
            }
            response.append(item)
            # print(scene['acquisition_date'].strftime('%Y-%m-%d'))
            # Write scene footprints to disk
            # fname = f"{scene['landsat_product_id']}.geojson"
            # print(fname)
            # with open(fname, "w") as f:
            #    json.dump(scene['spatial_coverage'].__geo_interface__, f)

        api.logout()

        return json.dumps(response)

    @app.route('/ee/download/<sceneID>')
    def download(sceneID):
        # Initialize a new API instance and get an access key
        api = EarthExplorer(app.config["USGS_USERNAME"],app.config["USGS_PASSWORD"])

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data")
      

        
        # Search for Landsat TM scenes
        filename = api.download(
            identifier = sceneID, 
            output_dir = DATA_FOLDER
        )

        ZIP_FILE = filename
        tar = tarfile.open(ZIP_FILE, "r:gz")
        tar.extractall(DATA_FOLDER + "/" + sceneID)
        tar.close()
        # os.remove(ZIP_FILE)

        # print (scenes)
        # print(f"{filename} scenes found.")

        api.logout()

        return filename

    @app.route('/ee/lst/<sceneID>')
    def lst(sceneID):
        # Initialize a new API instance and get an access key
        username = app.config["USGS_USERNAME"]
        password = app.config["USGS_PASSWORD"]
        api = API(username, password)

        # Collection 1
        product_id = api.get_display_id(sceneID, "landsat_8_c1")
        # assert product_id == "LT05_L1GS_173058_20111028_20161005_01_T2"
        # metadata = api.metadata(sceneID, "landsat_8_c1")

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data")
        # ZIP_FILE = DATA_FOLDER + "/" + product_id + ".tar.gz"
        METADATA_FILE = DATA_FOLDER + "/" + sceneID + "/" + product_id + "_MTL.txt"
        TEMP_FOLDER = DATA_FOLDER + "/landsat_lst/"

        LST_FOLDER = os.path.join(TEMP_FOLDER, sceneID)
        if not os.path.exists( LST_FOLDER ):
            _ = os.makedirs( LST_FOLDER )

        lst_retriever = LSTEstimator(METADATA_FILE, 'auto-ndvi-raw', LST_FOLDER)

        lst_retriever.get_lst_array()

        # print (scenes)
        # print(f"{filename} scenes found.")

        api.logout()

        return LST_FOLDER 

    # -u longlt_hcmut -P Longlt1870385 -v -p MOD11A2.006 -s MOLT -y 2020 -t h28v07 -o /home/netfpga/Workspace/remoteSensing/modisData/temp -b 265 -e 281
    @app.route('/modis/download/<year>/<startDate>/<endDate>')
    def modis_download(year, startDate = 0, endDate = 0):
        # Initialize a new API instance and get an access key
        username = app.config["EARTHDATA_USERNAME"]
        password = app.config["EARTHDATA_PASSWORD"]
        product = 'MOD11A2.006'
        platform = 'MOLT'
        tile = 'h28v07'
        URL= 'https://e4ftl01.cr.usgs.gov'
        proxy = None

        start = int(startDate)
        end = int(endDate)

        print(start)
        print(end)

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data/modis/")

        files = earthdata.get_modisfiles(username, password, platform, product, int(year), tile, proxy, doy_start = start, doy_end= end, base_url= URL, out_dir = DATA_FOLDER, verbose= False, ruff=False)
        
        print(f"{len(files)} found")

        # return username + "-" +  password + "-" + platform + "-" + product + "-" + tile + "-" + year + "-" + str(start) + "-" + str(end)
        return json.dumps(files)

    @app.route('/modis/lst/<sceneID>')
    def modis_lst(sceneID):
        # Initialize a new API instance and get an access key

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data/") 
        BASE_FOLDER = os.path.join(APP_ROOT, "data/modis/lst/") 
        SHP_FOLDER = DATA_FOLDER + "/HCM/"
        SHP_FILE = SHP_FOLDER + 'HCM.shp'
        TEMP_FOLDER = os.path.join(APP_ROOT, "data/modis_lst/")

        filename = sceneID + ".hdf"
        files = []
        files.append("%s%s" % (BASE_FOLDER, filename))

        print(files)


        DATA_FOLDER = os.path.join(TEMP_FOLDER, sceneID) 
        if not os.path.exists( DATA_FOLDER ):
            _ = os.makedirs( DATA_FOLDER )

        # SETUP OUTPUT DIRECTORIES
        converted_dir = os.path.join(DATA_FOLDER, 'converted')
        if not os.path.exists( converted_dir ):
            _ = os.makedirs( converted_dir )
        
        # make a tempdir to dump MOD11A2 bands
        temp_dir = os.path.join(DATA_FOLDER, 'TEMP')
        if not os.path.exists( temp_dir ):
            _ = os.makedirs( temp_dir )
        
        # make a directory to store the unwanted extra bands from the HDF files. (might be needed later)
        extra_bands_dir = os.path.join(DATA_FOLDER, 'raw_extra_bands_extracted')
        if not os.path.exists( extra_bands_dir ):
            _ = os.makedirs( extra_bands_dir )
        
        # make a directory for the mosaicked outputs
        mosaicked_dir = os.path.join( DATA_FOLDER, 'mosaicked' )
        if not os.path.exists( mosaicked_dir ):
            _ = os.makedirs( mosaicked_dir )
        
        # make a directory to store the warped to 32648 outputs
        warped_dir = os.path.join( DATA_FOLDER, 'warped' )
        if not os.path.exists( warped_dir ):
            _ = os.makedirs( warped_dir )
        
        # make a directory to store the final rescaled outputs
        rescaled_dir = os.path.join( DATA_FOLDER, 'rescaled' )
        if not os.path.exists( rescaled_dir ):
            _ = os.makedirs( rescaled_dir )
        
        # make a directory to store the final rescaled outputs
        out_dir = os.path.join( DATA_FOLDER, 'output' )
        if not os.path.exists( out_dir ):
            _ = os.makedirs( out_dir )

        # -----
        # PROCESS BAND EXTRACTION AND CONVERSION TO GTiff
        out = earthdata.run_all_extract_bands(files, temp_dir, ncpus=14 )
        print("Temp Dir is " + temp_dir)
        # move the desired bands to a new location
        _ = earthdata.move_desired_bands(temp_dir, converted_dir, keep_bands=['01','05'], ncpus=14 )
        
        # move the bands we want to the extra_bands_dir
        _ = earthdata.move_undesired_bands(temp_dir, extra_bands_dir, ncpus=14 )
        
        # ---
        # MOSAIC TILES: using gdal_merge.py
        files = glob.glob( os.path.join( converted_dir, '*.tif' ) )
        args = earthdata.make_mosaic_args(files, mosaicked_dir)
        mosaicked = earthdata.run_mosaic_tiles( args, ncpus=5 )
        
        # ---
        # WARP TO EPSG:32648 HCM
        files = glob.glob( os.path.join( DATA_FOLDER, 'mosaicked', '*.tif' ) )
        args = [(fn, os.path.join( warped_dir, os.path.basename(fn) )) for fn in files]
        warped = earthdata.run_warp_to_32648( args, ncpus=7 )
        
        # ---
        # RESCALE TO THE PROPER VALUES -- according to user-guide
        # 	https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
        files = glob.glob( os.path.join( DATA_FOLDER, 'warped', '*.tif' ) )
        rescaled = earthdata.run_rescale_values( files, ncpus=10 )

        out_files = [ os.path.basename(fn) for fn in rescaled ]

        for file in rescaled:
            outfile = os.path.join( out_dir, os.path.basename(file))
            command = "gdalwarp -t_srs EPSG:32648 -of GTiff -cutline " + SHP_FILE + " -cl HCM -crop_to_cutline -dstnodata 0.0 " + file + " " + outfile
            output = os.system(command)

        return json.dumps(out_files)


    @app.route('/modis/aod/download/<year>/<startDate>/<endDate>')
    def modis_aod_download(year, startDate = 0, endDate = 0):
        # Initialize a new API instance and get an access key
        username = app.config["EARTHDATA_USERNAME"]
        password = app.config["EARTHDATA_PASSWORD"]
        product = 'MCD19A2'
        tile = 'h28v07'
        URL= 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD19A2/' + year + "/" + startDate
        proxy = None

        start = int(startDate)
        end = int(endDate)

        print(start)
        print(end)

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data/modis_aod/")

        files = earthdata_aod.get_modisfiles(username, password, product, int(year), tile, proxy, doy_start = start, doy_end= end, base_url= URL, out_dir = DATA_FOLDER, verbose= False, ruff=False)
        
        print(f"{len(files)} found")

        # return username + "-" +  password + "-" + platform + "-" + product + "-" + tile + "-" + year + "-" + str(start) + "-" + str(end)
        return json.dumps(files)

    
    @app.route('/modis/aod/<sceneID>')
    def modis_aod(sceneID):
        # Initialize a new API instance and get an access key

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data/") 
        BASE_FOLDER = os.path.join(APP_ROOT, "data/modis/aod/") 
        SHP_FOLDER = DATA_FOLDER + "/HCM/"
        SHP_FILE = SHP_FOLDER + 'HCM.shp'
        TEMP_FOLDER = os.path.join(APP_ROOT, "data/modis_aod/") 


        filename = sceneID + ".hdf"
        files = []
        files.append("%s%s" % (BASE_FOLDER, filename))

        print(files)

        DATA_FOLDER = os.path.join(TEMP_FOLDER, sceneID) 
        if not os.path.exists( DATA_FOLDER ):
            _ = os.makedirs( DATA_FOLDER )

        # SETUP OUTPUT DIRECTORIES
        converted_dir = os.path.join(DATA_FOLDER, 'converted')
        if not os.path.exists( converted_dir ):
            _ = os.makedirs( converted_dir )
        
        # make a tempdir to dump MOD11A2 bands
        temp_dir = os.path.join(DATA_FOLDER, 'TEMP')
        if not os.path.exists( temp_dir ):
            _ = os.makedirs( temp_dir )
        
        # make a directory to store the unwanted extra bands from the HDF files. (might be needed later)
        extra_bands_dir = os.path.join(DATA_FOLDER, 'raw_extra_bands_extracted')
        if not os.path.exists( extra_bands_dir ):
            _ = os.makedirs( extra_bands_dir )
        
        # make a directory for the mosaicked outputs
        mosaicked_dir = os.path.join( DATA_FOLDER, 'mosaicked' )
        if not os.path.exists( mosaicked_dir ):
            _ = os.makedirs( mosaicked_dir )
        
        # make a directory to store the warped to 3338 outputs
        warped_dir = os.path.join( DATA_FOLDER, 'warped' )
        if not os.path.exists( warped_dir ):
            _ = os.makedirs( warped_dir )
        
        # make a directory to store the final rescaled outputs
        rescaled_dir = os.path.join( DATA_FOLDER, 'rescaled' )
        if not os.path.exists( rescaled_dir ):
            _ = os.makedirs( rescaled_dir )

        # make a directory to store the final rescaled outputs
        out_dir = os.path.join( DATA_FOLDER, 'output' )
        if not os.path.exists( out_dir ):
            _ = os.makedirs( out_dir )
        
        # data3, hdf,AOD_Uncertainty,Column_WV, optical_depth = earthdata_aod.processHDF(files[0])

        # modisData = data3.iloc[index,:]

        # modisData.to_csv("modisData.csv")
        # -----
        # PROCESS BAND EXTRACTION AND CONVERSION TO GTiff
        out = earthdata_aod.run_all_extract_bands(files, temp_dir, ncpus=14 )
        print("Temp Dir is " + temp_dir)

        # move the desired bands to a new location
        _ = earthdata_aod.move_desired_bands(temp_dir, converted_dir, keep_bands=['01','02'], ncpus=14 )
        
        # move the bands we want to the extra_bands_dir
        _ = earthdata_aod.move_undesired_bands(temp_dir, extra_bands_dir, ncpus=14 )
        
        # ---
        # MOSAIC TILES: using gdal_merge.py
        files = glob.glob( os.path.join( converted_dir, '*.tif' ) )
        args = earthdata_aod.make_mosaic_args(files, mosaicked_dir)
        mosaicked = earthdata_aod.run_mosaic_tiles( args, ncpus=5 )
        
        # ---
        # WARP TO EPSG:32648 HCM
        files = glob.glob( os.path.join( DATA_FOLDER, 'mosaicked', '*.tif' ) )
        args = [(fn, os.path.join( warped_dir, os.path.basename(fn) )) for fn in files]
        warped = earthdata_aod.run_warp_to_3338( args, ncpus=7 )
        
        # ---
        # RESCALE TO THE PROPER VALUES -- according to user-guide
        # 	https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
        files = glob.glob( os.path.join( DATA_FOLDER, 'warped', '*.tif' ) )
        rescaled = earthdata_aod.run_rescale_values( files, ncpus=10 )

        out_files = [ os.path.basename(fn) for fn in rescaled ]

        for file in rescaled:
             outfile = os.path.join( out_dir, os.path.basename(file))
             command = "gdalwarp -t_srs EPSG:32648 -of GTiff -cutline " + SHP_FILE + " -cl HCM -crop_to_cutline -dstnodata 0.0 " + file + " " + outfile
             output = os.system(command)

        return json.dumps(rescaled)

    @app.route('/ee/mosaic/<sceneList>')
    def mosaic(sceneList):
        # Initialize a new API instance and get an access key
        sceneArray = sceneList.split(',')
        print(sceneArray)

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data")
        # ZIP_FILE = DATA_FOLDER + "/" + product_id + ".tar.gz"
        TEMP_FOLDER = DATA_FOLDER + "/landsat_lst/"
        SHP_FOLDER = DATA_FOLDER + "/HCM/"
        SHP_FILE = SHP_FOLDER + 'HCM.shp'

        x = datetime.datetime.now()
        MOSAIC_FOLDER = os.path.join(TEMP_FOLDER, x.strftime("%Y%m%d"))
        if not os.path.exists( MOSAIC_FOLDER ):
            _ = os.makedirs( MOSAIC_FOLDER )
            for sceneID in sceneArray:
                LST_FOLDER = os.path.join(TEMP_FOLDER, sceneID)    
                LST_B10 =  sceneID + '_LST_B10'
                LST_B11 =  sceneID + '_LST_B11'
                # move the desired bands to a new location
                _ = copy_desired_tiff(LST_FOLDER, MOSAIC_FOLDER, keep_bands=[LST_B10,LST_B11], ncpus=14 )


        # Make a search criteria to select the DEM files
        search_B10 = "*_B10.tif"
        search_B11 = "*_B11.tif"

        q_b10 = os.path.join(MOSAIC_FOLDER, search_B10)
        q_b11 = os.path.join(MOSAIC_FOLDER, search_B11)

        # glob function can be used to list files from a directory with specific criteria
        B10List = glob.glob(q_b10)
        B11List = glob.glob(q_b11)

        # Files that were found:
        outVRT_B10 = MOSAIC_FOLDER + '/LC8_' + x.strftime("%Y%m%d") + '_LST_B10.vrt'
        outVRT_B11 = MOSAIC_FOLDER + '/LC8_' + x.strftime("%Y%m%d") + '_LST_B11.vrt'

        vrt_options = gdal.BuildVRTOptions(resampleAlg='average', addAlpha=False)

        b10_vrt = gdal.BuildVRT(outVRT_B10, B10List, options=vrt_options)
        b10_vrt = None

        b11_vrt = gdal.BuildVRT(outVRT_B11, B11List, options=vrt_options)
        b11_vrt = None

        outfile = MOSAIC_FOLDER + "/HCM_LST_B10.tif"
        command = "gdalwarp -s_srs EPSG:32648 -t_srs EPSG:32648 -of GTiff -cutline " + SHP_FILE + " -cl HCM -crop_to_cutline -dstnodata 0.0 " + outVRT_B10 + " " + outfile
        output = os.system(command)

        # print (scenes)
        # print(f"{filename} scenes found.")

        return send_file(outfile, mimetype='image/tiff', \
                                     as_attachment=True)

    @app.route('/modis/mosaic/<sceneList>')
    def mosaic_modis(sceneList):
        # Initialize a new API instance and get an access key
        sceneArray = sceneList.split(',')
        print(sceneArray)

        APP_ROOT = os.path.dirname(os.path.abspath(__file__))
        DATA_FOLDER = os.path.join(APP_ROOT, "data")
        # ZIP_FILE = DATA_FOLDER + "/" + product_id + ".tar.gz"
        TEMP_FOLDER = DATA_FOLDER + "/modis_lst/"
        SHP_FOLDER = DATA_FOLDER + "/HCM/"
        SHP_FILE = SHP_FOLDER + 'HCM.shp'

        x = datetime.datetime.now()
        MOSAIC_FOLDER = os.path.join(TEMP_FOLDER, x.strftime("%Y%m%d"))
        if not os.path.exists( MOSAIC_FOLDER ):
            _ = os.makedirs( MOSAIC_FOLDER )
            for sceneID in sceneArray:
                LST_FOLDER = os.path.join(TEMP_FOLDER, sceneID)    
                LST_B10 =  sceneID + '_LST_B10'
                LST_B11 =  sceneID + '_LST_B11'
                # move the desired bands to a new location
                _ = copy_desired_tiff(LST_FOLDER, MOSAIC_FOLDER, keep_bands=[LST_B10,LST_B11], ncpus=14 )


        # Make a search criteria to select the DEM files
        search_B10 = "*_B10.tif"
        search_B11 = "*_B11.tif"

        q_b10 = os.path.join(MOSAIC_FOLDER, search_B10)
        q_b11 = os.path.join(MOSAIC_FOLDER, search_B11)

        # glob function can be used to list files from a directory with specific criteria
        B10List = glob.glob(q_b10)
        B11List = glob.glob(q_b11)

        # Files that were found:
        outVRT_B10 = MOSAIC_FOLDER + '/LC8_' + x.strftime("%Y%m%d") + '_LST_B10.vrt'
        outVRT_B11 = MOSAIC_FOLDER + '/LC8_' + x.strftime("%Y%m%d") + '_LST_B11.vrt'

        vrt_options = gdal.BuildVRTOptions(resampleAlg='average', addAlpha=False)

        b10_vrt = gdal.BuildVRT(outVRT_B10, B10List, options=vrt_options)
        b10_vrt = None

        b11_vrt = gdal.BuildVRT(outVRT_B11, B11List, options=vrt_options)
        b11_vrt = None

        outfile = MOSAIC_FOLDER + "/HCM_LST_B10.tif"
        command = "gdalwarp -s_srs EPSG:32648 -t_srs EPSG:32648 -of GTiff -cutline " + SHP_FILE + " -cl HCM -crop_to_cutline -dstnodata 0.0 " + outVRT_B10 + " " + outfile
        output = os.system(command)

        # print (scenes)
        # print(f"{filename} scenes found.")

        return send_file(outfile, mimetype='image/tiff', \
                                     as_attachment=True)

    return app