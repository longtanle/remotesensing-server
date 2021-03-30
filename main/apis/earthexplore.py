from flask_restplus import Namespace, Resource, fields

api = Namespace("EarthExplorer", description="EarthExplorer related APIs")
from flask import request
from flask_jwt_extended import jwt_required

"""Tests for ee_api module."""
import json

import os
from datetime import datetime
from shapely.geometry import Polygon
from core import ee_apis, ee_errors, ee_utils

from core.ee_earthexplorer import EarthExplorer
from core.ee_errors import EarthExplorerError

from main.services.earthexplore_service import EarthExplorerService

earthexplorer_model = api.model(
    "EarthExplorerModel",
    {
        "landsat_product_id": fields.String(description="Landsat 8 scence name", required=True),
        "acquisition_date": fields.String(description="Landsat 8 acquisition name", required=True),

        # "author": fields.String(description="Name of the author", required=True),
        #"genres": fields.String(description="Type of book", required=True),
        # "year": fields.String(description="year of publication", required=True),
    },
)

