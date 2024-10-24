# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 21:07:04 2024

@author: gabriel.ferraz
"""

# from planet import Session
import os
from os.path import isfile, join
import pandas as pd
import geopandas as gpd
import requests
from shapely.geometry import mapping
from shapely.geometry import Polygon
from shapely.geometry import shape
import rasterio
from rasterio import mask
from io import BytesIO
import io
import numpy as np
# from rasterio.mask import mask
import time
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
from itertools import islice
import datetime
import json
import logging
import boto3
# gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

os.chdir(r'C:\Projetos\bioflore\gbs-madagascar')
         
load_dotenv(find_dotenv())

#API Key Planet
API_KEY = os.getenv('PL_API_KEY')

# Specify your access key and secret key
ACCESS_KEY = os.getenv('ACCESS_KEY')
SECRET_KEY = os.getenv('SECRET_KEY')

global directory
directory = r"C:\Projetos\bioflore\gbs-madagascar"
now = datetime.datetime.today().strftime('%Y-%m-%dT%H%M%S.%fZ')

logging.basicConfig(
     filename=directory+f"/log_{now}.txt",
     filemode='a',
     format='%(asctime)s\t%(levelname)s\t%(message)s',
     level=logging.INFO
 )

# Define our AOI (Area of Interest)
# Read the shapefile
gdf = gpd.read_file(r"C:\Projetos\bioflore\gbs-madagascar\geo\Andohahela_parcel_I_&_II\Andohahela_parcel_I_&_II.shp")
# gdf2 = gpd.read_file(directory + "/geo/Contract174_Polygons.gpkg")
gdf = gdf.head(1)


# rdf = gpd.GeoDataFrame( pd.concat( [gdf1, gdf2], ignore_index=True) )
gdf.crs = "EPSG:4326"

logging.info('Initializing process...')

def remove_z_gen_aoi(geometry: gpd.geodataframe.GeoDataFrame) -> dict:
    try:
        # geometry = shape['geometry'][0]
        exterior = geometry.exterior.coords[:-1]  # Drop the z value for exterior
        flattened_exterior = [(x, y) for x, y, z in exterior]  # Keep only (x, y)
        
        flattened_interiors = []
        for interior in geometry.interiors:
            interior_coords = interior.coords[:-1]  # Drop the z value for interior
            flattened_interior = [(x, y) for x, y, z in interior_coords]
            flattened_interiors.append(flattened_interior)
        
        # Create a new 2D Polygon
        flattened_polygon = Polygon(flattened_exterior, flattened_interiors)
        # Convert the GeoJSON geometry to a Python dict
        aoi = mapping(flattened_polygon)
        return aoi
    except Exception as e:
       logging.error(f'Treat Polygon: {e}')
       
def gen_aoi(geometry: gpd.geodataframe.GeoDataFrame) -> dict:
    try:
        # geometry = gdf['geometry'][0]
        exterior = geometry.exterior.coords[:-1]  # Drop the z value for exterior
        flattened_exterior = [(x, y) for x, y in exterior]  # Keep only (x, y)
        
        flattened_interiors = []
        for interior in geometry.interiors:
            interior_coords = interior.coords[:-1]  # Drop the z value for interior
            flattened_interior = [(x, y) for x, y, z in interior_coords]
            flattened_interiors.append(flattened_interior)
        
        # Create a new 2D Polygon
        flattened_polygon = Polygon(flattened_exterior, flattened_interiors)
        # Convert the GeoJSON geometry to a Python dict
        aoi = mapping(flattened_polygon)
        return aoi
    except Exception as e:
       logging.error(f'Treat Polygon: {e}')

def request_planet_ids(aoi : dict, gte : str, lte: str,\
                       cloud: float, item_types: list, ) -> pd.DataFrame:
    try:
        #request Planet API
        # filter for items the overlap with our chosen geometry
        # aoi = aoi
        # gte ='2019-01-01T00:00:00.000Z'
        # lte = '2024-07-01T00:00:00.000Z'
        # cloud = 0.01
        # item_types = ['PSScene']
        
        geometry_filter = {
          "type": "GeometryFilter",
          "field_name": "geometry",
          "config": aoi
        }
        
        # filter images acquired in a certain date range
        date_range_filter = {
          "type": "DateRangeFilter",
          "field_name": "acquired",
          "config": {
            "gte": gte,
            "lte": lte
          }
        }
        
        # filter any images which are more than x% clouds
        cloud_cover_filter = {
          "type": "RangeFilter",
          "field_name": "cloud_cover",
          "config": {
            "lte": cloud
          }
        }
        
        permission_filter= {
               "type":"PermissionFilter",
               "config":[
                  "assets:download"
               ]
            }
        
        not_filter = {
       "type":"NotFilter",
       "config":{
             "type":"StringInFilter",
             "field_name":"quality_category",
             "config":[
                "test"
             ]
         }
    }
        asset_filter =   {
                    "type": "AssetFilter",
                    "config": [
                        "ortho_analytic_4b_sr"
                    ]
                }
        # create a filter that combines our geo and date filters
        # could also use an "OrFilter"
        redding_reservoir = {
          "type": "AndFilter",
          "config": [geometry_filter, date_range_filter, cloud_cover_filter,
                     permission_filter, not_filter, asset_filter]
        }
        
        
        # Stats API request object
        search_endpoint_request = {
          "item_types": item_types,
          "filter": redding_reservoir,
          
        }
        
        session = requests.Session()
        session.auth = (API_KEY, '')
        
        # fire off the POST request
        iterate_pages = []
        
        result = \
          session.post(
             # 'https://api.planet.com/data/v1/stats',
              'https://api.planet.com/data/v1/quick-search?_sort=acquired asc&_page_size=250',
            # auth=HTTPBasicAuth(API_KEY, ''),
            json=search_endpoint_request)
        
        data = result.json()
        
        # grid_cell = data['features'][0]['properties']['grid_cell']
        ids = [(f['id'], f['properties']['acquired'], f['properties']['cloud_cover'],f['geometry'])\
               for f in data['features']]
            
        iterate_pages.extend(ids)
        
        while '_next' in data['_links']:
            try:
                result = \
                  session.get(
                     # 'https://api.planet.com/data/v1/stats',
                      data['_links']['_next'])
                    # auth=HTTPBasicAuth(API_KEY, ''),
                    # json=search_endpoint_request)
                data = result.json()
                
                # grid_cell = data['features'][0]['properties']['grid_cell']
                ids = [(f['id'], f['properties']['acquired'], f['properties']['cloud_cover'],f['geometry'])\
                       for f in data['features']]
                iterate_pages.extend(ids)
            except:
                break
            
        # unique_tuples = list(set(iterate_pages))
        df = pd.DataFrame(iterate_pages, columns=['id', 'date', 
                                                  'cloud', 'geom'])
    
        df['date'] = pd.to_datetime(df['date'])
        df.sort_values(by='date', inplace=True)
        # df['date'] = df['date'].dt.date
    
        # Drop duplicate rows based on the modified 'date' column
        # df = df.drop_duplicates(subset=['date'], keep='first')
        
        return df
    except Exception as e:
        logging.error(f'Quick Search: {e}')  
    
def check_intersection(row):
    return shape(row['geom']).contains_properly(limite)

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_pairs(data):
    it = iter(data)
    for i in range(0, len(data), 2):
        yield {k:data[k] for k in islice(it, 2)}
            
def request_download_urls_planet(ids_list, item_type, asset_type,
                                 aoi, order_name):
    try:
        # ids_list = chunked_ids[0]
        
        session = requests.Session()
        session.auth = (API_KEY, '')
        
        call = {
           "name": order_name,
           "source_type": "scenes",
           "products":[
              {
                 "item_ids": ids_list,
                 "item_type":item_type,
                 "product_bundle":asset_type
              }
           ],
       "tools": [
           {
               "clip": {
                   "aoi": aoi
                   }
               }
           ]
        }
      
        req = session.post( 'https://api.planet.com/compute/ops/orders/v2', 
                          json=call)
        
        data_req = req.json()
        order_id = data_req['id']
        # order_id = '643ad5dc-3c4e-45e5-9d34-0f380efe4f80'
        
        # request an item
        url_req =  f"https://api.planet.com/compute/ops/orders/v2/{order_id}"
        
        item = \
          session.get(url_req
           )
        
        data_dwnld = item.json()
        
        while 'results' not in data_dwnld['_links']:
            logging.info('Waiting URLs to be available')
            time.sleep(120)
            # print('still waiting...')
            item = \
              session.get(url_req
               )
        
            data_dwnld = item.json()
        
        dwnld_dict = {}
            
        for i in data_dwnld['_links']['results']:
            if i['name'].endswith('.tif'):
                dwnld_dict[i['name']] = i['location']
        
        return dwnld_dict
    except Exception as e:
        logging.error(f'Download Request: {e}')  
     
def download_from_url(url_name, talhao):
    try:
        # url_name = pairs[0]
        urls = list(url_name.values())
        date = list(url_name.keys())[0].split('/')[-1][:15]
        img_name =  f'{talhao}_{date}.tif'
        qa_name = f'{talhao}_{date}_QA.tif'
        
        response1 = requests.get(urls[1])
        img_data1 = response1.content
        with open(directory +'/PSData/'+ img_name, 'wb') as handler:
            handler.write(img_data1)
        # s3.Bucket(BUCKET_NAME).upload_file(directory +'/'+ img_name,
        #                               s3_folder + img_name)
        # src1 = rasterio.open(BytesIO(response1.content))
        # sr_data = src1.read()

        response2 = requests.get(urls[0])
        img_data2 = response2.content
        with open(directory +'/PSData/'+ qa_name, 'wb') as handler:
            handler.write(img_data2)
        # s3.Bucket(BUCKET_NAME).upload_file(directory +'/'+ qa_name,
        #                               s3_folder +  qa_name)
        # src2 = rasterio.open(BytesIO(response2.content))
        # qa_data = src2.read(1)
        
        # mask_array = np.where(qa_data == 0, 0, 1)
        
        # Apply the mask to the surface reflectance raster
        # masked_sr_data = sr_data * mask_array
        
        # nir_band = masked_sr_data[3]
        # red_band = masked_sr_data[2]
        # ndvi = (nir_band - red_band) / (nir_band + red_band)
        # ndvi_mean = np.nanmean(ndvi, dtype=np.float64)
        
        # result = (date, ndvi_mean, talhao)
        
        # os.remove(directory +'/'+ img_name)
        # os.remove(directory +'/'+ qa_name)
        
        return
    except Exception as e:
        logging.error(f'Get NDVI: {e}')  
        
for limite, talhao in gdf[['geometry', 'SITE_NAME1',
                            ]].itertuples(index=False):

    # comment later
    # limite = gdf['geometry'][0]
    # talhao = 'extension'
        # continue
    plantig_date = '2022-01-01'
    
    item_type = 'PSScene'
    asset_type = "analytic_sr_udm2"
    
    logging.info(f'***********Starting parcel {talhao}***********')
    
    if not limite.is_valid:
        logging.info(f'Polygon {talhao} not valid. Correction applied')
        limite = limite.buffer(0)
        
    aoi = gen_aoi(limite)
    
    sd = datetime.datetime.strptime(plantig_date, "%Y-%m-%d").strftime(
        '%Y-%m-%dT%H:%M:%S.%fZ')
    ed = datetime.datetime.today().strftime(
        '%Y-%m-%dT%H:%M:%S.%fZ')
    
    ids = request_planet_ids(aoi, sd,\
                       ed, 0.1, [item_type])
    
    # if ids is None:
    #     continue
    # Create a new column 'fully_intersects' with boolean values
    ids['fully_intersects'] = ids.apply(check_intersection, axis=1)
    
        
    # Filter out rows that do not fully intersect with the polygon
    filtered_ids = ids[ids['fully_intersects']]
    idx_min = filtered_ids.groupby(filtered_ids['date'].dt.to_period('M'))['cloud'].idxmin()
    r = filtered_ids.loc[idx_min]
    
    # len(r)
    
    chunked_ids = list(chunks(r['id'].tolist(), 500))
       
    dwnld_dict = {}
    for i in chunked_ids:
        # i = chunked_ids[0]
        with open(directory + "/ids_list.txt", 'w') as file:
            file.write(str(i))
        now = datetime.datetime.today().strftime('%Y-%m-%dT%H%M%S.%fZ')
        order_name = f'talhao_{talhao}_{now}'
        dwnld_dict.update(request_download_urls_planet(i, item_type, asset_type,
                                                       aoi, order_name))
    with open(directory + "/donwload_dict.txt", 'r') as file:
        data = file.read().replace('\n', '')
    
    pairs = list(get_pairs(dwnld_dict))
    
    for pair in pairs:
        download_from_url(pair, talhao)
            
    
    logging.info('Process Complete!')    
    
  