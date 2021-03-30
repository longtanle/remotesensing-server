import datetime
import logging as log
import os

import coloredlogs
from flask import Flask
from flask_bcrypt import Bcrypt
from flask_cors import CORS
from flask_jwt_extended import (
    JWTManager,
    create_access_token,
    get_jwt_identity,
    jwt_required,
)

import subprocess, glob, shutil

import json
from core.ee_apis import API
import tempfile
from core.ee_earthexplorer import EarthExplorer, EarthExplorerError
import core.earthdata as earthdata
import core.earthdata_aod as earthdata_aod
from core.LSTEstimator import LSTEstimator

import zipfile
import tarfile



# local imports
from config import app_config

from .db import MongoDB


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
        return "Welcome to Remote Sensing Server!"

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
        os.remove(ZIP_FILE)

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
        ZIP_FILE = DATA_FOLDER + "/" + product_id + ".tar.gz"
        METADATA_FILE = DATA_FOLDER + "/" + sceneID + "/" + product_id + "_MTL.txt"
        TEMP_FOLDER = os.path.join(APP_ROOT, "data/temp")

        lst_retriever = LSTEstimator(METADATA_FILE, 'auto-ndvi-raw',TEMP_FOLDER)

        lst_retriever.get_lst_array()

        # print (scenes)
        # print(f"{filename} scenes found.")

        api.logout()

        return ZIP_FILE 

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
        BASE_FOLDER = os.path.join(APP_ROOT, "data/modis/") 

        filename = sceneID + ".hdf"
        files = []
        files.append("%s%s" % (BASE_FOLDER, filename))

        print(files)


        DATA_FOLDER = os.path.join(BASE_FOLDER, sceneID) 
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
        # WARP TO EPSG:3338 AK Albers 
        files = glob.glob( os.path.join( DATA_FOLDER, 'mosaicked', '*.tif' ) )
        args = [(fn, os.path.join( warped_dir, os.path.basename(fn) )) for fn in files]
        warped = earthdata.run_warp_to_3338( args, ncpus=7 )
        
        # ---
        # RESCALE TO THE PROPER VALUES -- according to user-guide
        # 	https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
        files = glob.glob( os.path.join( DATA_FOLDER, 'warped', '*.tif' ) )
        rescaled = earthdata.run_rescale_values( files, ncpus=10 )

        return json.dumps(files)


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
        BASE_FOLDER = os.path.join(APP_ROOT, "data/modis_aod/") 


        filename = sceneID + ".hdf"
        files = []
        files.append("%s%s" % (BASE_FOLDER, filename))

        print(files)

        DATA_FOLDER = os.path.join(BASE_FOLDER, sceneID) 
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
        mosaicked = earthdata.run_mosaic_tiles( args, ncpus=5 )
        
        # ---
        # WARP TO EPSG:3338 AK Albers 
        files = glob.glob( os.path.join( DATA_FOLDER, 'mosaicked', '*.tif' ) )
        args = [(fn, os.path.join( warped_dir, os.path.basename(fn) )) for fn in files]
        warped = earthdata_aod.run_warp_to_3338( args, ncpus=7 )
        
        # ---
        # RESCALE TO THE PROPER VALUES -- according to user-guide
        # 	https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
        files = glob.glob( os.path.join( DATA_FOLDER, 'warped', '*.tif' ) )
        rescaled = earthdata_aod.run_rescale_values( files, ncpus=10 )

        return json.dumps(files)

    return app