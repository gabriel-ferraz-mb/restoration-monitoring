# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 19:31:37 2024

@author: Gabriel
"""

import requests
import pandas
import json
import os

os.chdir(r'C:\Projetos\bioflore\gbs-uganda\geo\gfw')

# =============================================================================
# CREATE USER
# =============================================================================
params ={
"name": "Gabriel Ferraz",
"email": "gabriel.rezendeferraz@gmail.com"
}

headers = {
        'Content-Type': 'application/json'
        }

response = requests.post('https://data-api.globalforestwatch.org/auth/sign-up',
                         headers=headers, json=params)
print(response.json())

with open('sign-up.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)

# =============================================================================
# GET TOKEN
# =============================================================================
params ={
"username": "gabriel.rezendeferraz@gmail.com",
"password": "Capio236!"
}

headers = {
        'Content-Type': 'application/x-www-form-urlencoded'
        }

response = requests.post('https://data-api.globalforestwatch.org/auth/token',
                         headers=headers, data=params)
print(response.json())

with open('token.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)

token = response.json()['data']['access_token']

# =============================================================================
# GET API KEY
# =============================================================================
params ={
 "alias": "api-key-bioflore",
 "email": "gabriel.rezendeferraz@gmail.com.br",
 "organization": "Bioflore",
 "domains": []
}

headers = {
        'Authorization' : f'Bearer {token}',
        'Content-Type': 'application/json'
        }

response = requests.post('https://data-api.globalforestwatch.org/auth/apikey',
                         headers=headers, json=params)
print(response.json())

with open('api-key.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)

user_id = response.json()['data']['user_id']
api_key = response.json()['data']['api_key']

# =============================================================================
# GET DATASETS
# =============================================================================

headers = {
        'x-api-key': api_key
        }

response = requests.get('https://data-api.globalforestwatch.org/datasets',
                         headers=headers)
print(response.json())

with open('datasets.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)

# =============================================================================
# GET GFW INTEGRATED ALERTS METADATA
# =============================================================================

headers = {
        'x-api-key': api_key
        }


response = requests.get('https://data-api.globalforestwatch.org/dataset/gfw_integrated_alerts',
                         headers=headers)

print(response.json())

with open('gfw_integrated_alerts.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)
    
# =============================================================================
# GET GFW INTEGRATED ALERTS METADATA FIELDS
# =============================================================================

headers = {
        'x-api-key': api_key
        }


response = requests.get('https://data-api.globalforestwatch.org/dataset/gfw_integrated_alerts/latest/fields',
                         headers=headers)

print(response.json())

with open('gfw_integrated_alerts_fields.json', 'w', encoding='utf-8') as f:
    json.dump(response.json(), f, ensure_ascii=False, indent=4)

# =============================================================================
# GET DATASET
# =============================================================================
import geopandas as gpd
from shapely.geometry import box
import ast
import pandas as pd

def get_extent_as_geodf(gdf):
    # Get the total bounds of the GeoDataFrame
    bounds = gdf.total_bounds  # Returns an array: [minx, miny, maxx, maxy]

    # Create a bounding box (rectangle) using shapely
    bbox = box(bounds[0], bounds[1], bounds[2], bounds[3])

    # Create a new GeoDataFrame containing just this bounding box
    extent_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs=gdf.crs)

    return extent_gdf

api_key = '00624f70-1ec7-4631-8b4b-4d4f75b2bbe3'

polygons = gpd.read_file("../shp/restoration_sites.shp")
extent = get_extent_as_geodf(polygons)
j = ast.literal_eval(extent.to_json())

geometry = j['features'][0]['geometry']
q = 'SELECT longitude, latitude, gfw_integrated_alerts__date,gfw_integrated_alerts__confidence from results'

params= {
        'geometry':geometry,
        'sql' : q
        }

headers = {
        'x-api-key': api_key,
        'Content-Type': 'application/json' 
        }


response = requests.post('https://data-api.globalforestwatch.org/dataset/gfw_integrated_alerts/latest/query',
                         headers=headers, json = params)


results = pd.DataFrame(response.json()['data'])

gdf = gpd.GeoDataFrame(
   results, geometry=gpd.points_from_xy(results.longitude, results.latitude), crs="EPSG:4326"
)

gdf['gfw_integrated_alerts__date'] = pd.to_datetime(gdf['gfw_integrated_alerts__date'], format="%Y-%m-%d")
gdf['year'] = pd.DatetimeIndex(gdf['gfw_integrated_alerts__date'] ).year.astype(str)

import matplotlib.pyplot as plt

# =============================================================================
# PLOT MAP
# =============================================================================

ig, ax = plt.subplots(figsize=(10, 10))

# Plot the polygons
polygons.boundary.plot(ax=ax, linewidth=1, color='black')

# Plot the points with color based on the 'year' column
gdf.plot(ax=ax, column='year', cmap='viridis', legend=True, 
         markersize=10, alpha=1)

# Add basemap if needed (optional, requires contextily)
# import contextily as ctx
# ctx.add_basemap(ax, crs=gdf.crs.to_string())

# Add title and labels, adjust plot limits
ax.set_title("Forest Disturbance by Year")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.show()

# =============================================================================
# PLOT CHART
# =============================================================================

year_counts = gdf['year'].value_counts().sort_index()
year_counts = year_counts*0.01

# Create the bar graph
fig, ax = plt.subplots(figsize=(10, 6))

year_counts.plot(kind='bar', ax=ax, color='skyblue')

# Set labels and title
ax.set_title("Hectares deforestated per Year")
ax.set_xlabel("Year")
ax.set_ylabel("Hectares")
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Optionally rotate x-axis labels if there are many years
plt.xticks(rotation=45)

plt.tight_layout()  # To adjust layout so everything fits
plt.show()