# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 20:45:38 2024

@author: Gabriel
"""
import eemont, geemap
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import seaborn as sns
from matplotlib import pyplot as plt
from datetime import date
from time import strptime
import geopandas as gpd
from shapely.geometry import mapping
from shapely.geometry import box
import warnings
warnings.filterwarnings('ignore')
import datetime
from datetime import time, timedelta, date
from datetime import datetime
from dateutil.relativedelta import relativedelta
import json
from google.oauth2 import service_account

def get_earth_engine_connection(key: str):
    credentials = service_account.Credentials.from_service_account_file(key)
    scoped_credentials = credentials.with_scopes(
        ['https://www.googleapis.com/auth/cloud-platform'])
    return scoped_credentials

def remove_z(geometry):
#     if geometry.is_empty:
#         return geometry  # Handle empty geometries if necessary
    # Convert geometry to 2D
    def to_2d(geom):
        if isinstance(geom, Polygon):
            # Extract 2D coordinates from each part (exterior and interiors)
            exterior_2d = [(x, y) for x, y, *_ in geom.exterior.coords]
            interiors_2d = [
                [(x, y) for x, y, *_ in interior.coords]
                for interior in geom.interiors
            ]
            return Polygon(exterior_2d, interiors_2d)
        elif isinstance(geom, MultiPolygon):
            # Iterate over individual Polygons within the MultiPolygon
            return MultiPolygon([to_2d(p) for p in geom.geoms])
        else:
            return geom

    return to_2d(geometry)

def geometry_to_ee(aoi):
    if aoi.geometry.unary_union.geom_type != 'Polygon':
        # mount the bounding box of the geometries
        clip_geometry = mapping(aoi.unary_union.convex_hull.buffer(0))
    else:
        clip_geometry = mapping(aoi.unary_union)
    # get coordinates
    coordinates = []
    for polygon in clip_geometry['coordinates']:
        poly_coords = []
        for coord in polygon:
            poly_coords.append(list(coord))
        coordinates.append(poly_coords)
    polygon = ee.Geometry.Polygon(coordinates)
    return polygon

def get_extent_as_geodf(gdf):
    # Get the total bounds of the GeoDataFrame
    bounds = gdf.total_bounds  # Returns an array: [minx, miny, maxx, maxy]

    # Create a bounding box (rectangle) using shapely
    bbox = box(bounds[0], bounds[1], bounds[2], bounds[3])

    # Create a new GeoDataFrame containing just this bounding box
    extent_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs=gdf.crs)

    return extent_gdf

import ee
KEY = r"C:\Projetos\bioflore\restoration_monitoring\bioflore-ee.json"
session = get_earth_engine_connection(KEY)
ee.Initialize(session)


def get_temperature_dataset(start_date, end_date, CAR, lat_long):
    point = ee.Geometry.Point(lat_long.getInfo()['coordinates'])#.buffer(buf)
    temp = ee.ImageCollection('NOAA/CFSR').select("Temperature_surface").filterBounds(CAR).filterDate(start_date, end_date)
    ts_temp = temp.getTimeSeriesByRegion(reducer = [ee.Reducer.mean()],
                                  geometry = point,
                                  bands = 'Temperature_surface',
                                  scale = 55660)
    tsPandas = geemap.ee_to_df(ts_temp)

    tsPandas_clean = tsPandas.loc[tsPandas['Temperature_surface'] > 0 ]
    return tsPandas_clean

def show_temperature_data(CAR, lat_long, years_past):
    
    # CAR = polygon
    # lat_long = lat_long
    # years_past = 2
    
    ts_tem_list = []
    for i in range(1,years_past):
        try:
            ts_tem = get_temperature_dataset(f'{str(2024-i)}-01-01', f'{str(2024-i)}-12-31', CAR, lat_long)
            ts_tem_list.append(ts_tem)
        except:
            print(f"year: {str(2024-i)} failed")
            continue
    
    result_temp = pd.concat(ts_tem_list)

    result_temp["Temperatura_celsius"] = result_temp["Temperature_surface"] - 273

    result_temp['Ano'] = pd.DatetimeIndex(result_temp['date']).year
    result_temp['month'] = pd.DatetimeIndex(result_temp['date']).month
    result_temp['id'] = result_temp['month']

    result_temp.loc[result_temp.month == 1, "month"] = "Jan"
    result_temp.loc[result_temp.month == 2, "month"] = "Feb"
    result_temp.loc[result_temp.month == 3, "month"] = "Mar"
    result_temp.loc[result_temp.month == 4, "month"] = "Apr"
    result_temp.loc[result_temp.month == 5, "month"] = "May"
    result_temp.loc[result_temp.month == 6, "month"] = "Jun"
    result_temp.loc[result_temp.month == 7, "month"] = "Jul"
    result_temp.loc[result_temp.month == 8, "month"] = "Aug"
    result_temp.loc[result_temp.month == 9, "month"] = "Sep"
    result_temp.loc[result_temp.month == 10, "month"] = "Oct"
    result_temp.loc[result_temp.month == 11, "month"] = "Nov"
    result_temp.loc[result_temp.month == 12, "month"] = "Dec"

    result_temp_ano = result_temp.groupby(['Ano', 'month', 'id'])["Temperatura_celsius"].mean()
    result_temp_ano = result_temp_ano.reset_index()
    result_temp_ano.sort_values(by='Ano', ascending=True, inplace=True)
    result_temp_ano_min = result_temp_ano.groupby(['month','id'])["Temperatura_celsius"].min()
    result_temp_ano_max = result_temp_ano.groupby(['month','id'])["Temperatura_celsius"].max()
    result_temp_ano_mean = result_temp_ano.groupby(['month','id'])["Temperatura_celsius"].mean()
    df_temp = pd.concat([result_temp_ano_mean, result_temp_ano_min,result_temp_ano_max],
                        axis=1).reset_index()
    df_temp.columns = ['mes','id','T_MEAN', 'T_MIN', 'T_MAX']
    
    # df_temp['id'] = df_temp.apply(lambda row: strptime(row.mes,'%b').tm_mon , axis = 1)
    
    df_temp_ano_mean = result_temp.groupby('Ano')["Temperatura_celsius"].mean()
    df_temp_ano_max = result_temp_ano.groupby(['Ano'])["Temperatura_celsius"].max()
    df_temp_ano_min = result_temp_ano.groupby(['Ano'])["Temperatura_celsius"].min()
    df_temp_ano = pd.concat([df_temp_ano_mean, df_temp_ano_min,df_temp_ano_max],
                        axis=1).reset_index()
    
    df_temp_ano.columns = ['ano','T_MEAN', 'T_MIN', 'T_MAX']

    return df_temp, df_temp_ano 

def get_precipitation_dataset(start_date, end_date, CAR, lat_long):
    
    # start_date = '2022-07-01'
    # end_date = '2023-01-01'
    # CAR = polygon
    # lat_long = lat_long
    
    ppt = ee.ImageCollection("JAXA/GPM_L3/GSMaP/v6/operational").select(
        "hourlyPrecipRate").filterBounds(CAR).filterDate(start_date, end_date)
    point = ee.Geometry.Point(lat_long.getInfo()['coordinates'])#.buffer(20000)

    ts_ppt = ppt.getTimeSeriesByRegion(reducer = [ee.Reducer.sum()],
                                  geometry = point,
                                  bands = 'hourlyPrecipRate',
                                  scale = 100)

    tsPandas = geemap.ee_to_df(ts_ppt)
    return tsPandas

def show_precipitation_data( CAR, lat_long, years_past):
    
    # CAR = polygon
    # lat_long = lat_long
    # years_past = 2
    
    ts_ppt_list = []
    for i in range(1,years_past):
        ts_ppt = get_precipitation_dataset(f'{str(2024-i)}-01-01', f'{str(2024-i)}-12-31', CAR, lat_long)
        ts_ppt_list.append(ts_ppt)
        
    tsPandas= pd.concat(ts_ppt_list)
    
    tsPandas_clean = tsPandas.loc[tsPandas['hourlyPrecipRate'] >= 0 ]
    tsPandas_clean['date'] = pd.to_datetime(tsPandas_clean['date'])
    tsPandas_clean['Ano'] = pd.DatetimeIndex(tsPandas_clean['date']).year
    tsPandas_clean['month'] = pd.DatetimeIndex(tsPandas_clean['date']).month
    
    tsPandas_clean = tsPandas_clean.reindex(columns=['Ano','month', 'reducer','date','hourlyPrecipRate'])
    tsPandas_clean = tsPandas_clean[['Ano', 'month','hourlyPrecipRate']]
    tsPandas_clean.sort_values(by='month', ascending=True, inplace= True)
    
    result_ppt_sum = tsPandas_clean.groupby(by=[tsPandas_clean.Ano, tsPandas_clean.month]).sum().reset_index()
    result_ppt_avg = result_ppt_sum.groupby(by=[result_ppt_sum.month])['hourlyPrecipRate'].mean().rename("avg")
    result_ppt_min = result_ppt_sum.groupby(by=[result_ppt_sum.month])['hourlyPrecipRate'].min().rename("min")
    result_ppt_max = result_ppt_sum.groupby(by=[result_ppt_sum.month])['hourlyPrecipRate'].max().rename("max")
    result_ppt_cumsum = result_ppt_avg.cumsum().rename("sum_avg")
    result_ppt_cumsum_min = result_ppt_min.cumsum().rename("sum_min")
    result_ppt_cumsum_max = result_ppt_max.cumsum().rename("sum_max")
    ppt_df = pd.concat([result_ppt_avg, result_ppt_min, result_ppt_max, result_ppt_cumsum, result_ppt_cumsum_min, result_ppt_cumsum_max], axis=1)
    
    ppt_df = ppt_df.reset_index()
    
    ppt_df.loc[ppt_df.month == 1, "month"] = "Jan"
    ppt_df.loc[ppt_df.month == 2, "month"] = "Feb"
    ppt_df.loc[ppt_df.month == 3, "month"] = "Mar"
    ppt_df.loc[ppt_df.month == 4, "month"] = "Apr"
    ppt_df.loc[ppt_df.month == 5, "month"] = "May"
    ppt_df.loc[ppt_df.month == 6, "month"] = "Jun"
    ppt_df.loc[ppt_df.month == 7, "month"] = "Jul"
    ppt_df.loc[ppt_df.month == 8, "month"] = "Aug"
    ppt_df.loc[ppt_df.month == 9, "month"] = "Sep"
    ppt_df.loc[ppt_df.month == 10, "month"] = "Oct"
    ppt_df.loc[ppt_df.month == 11, "month"] = "Nov"
    ppt_df.loc[ppt_df.month == 12, "month"] = "Dec"
    
    ppt_df_year = tsPandas_clean.groupby(by=tsPandas_clean.Ano)['hourlyPrecipRate'].sum().reset_index()
    
    return ppt_df, ppt_df_year

aoi = gpd.read_file(r"C:\Projetos\bioflore\gbs-uganda\geo\shp\restoration_sites.shp")
extent = get_extent_as_geodf(aoi)
aoi_2d = remove_z(extent.to_crs(4326))

polygon = geometry_to_ee(aoi_2d)

lat_long = polygon.centroid()

ppt_df, ppt_df_year = show_precipitation_data(polygon, lat_long, 2)
temp_df, temp_df_year = show_temperature_data(polygon, lat_long, 2)

        
def multi_plot_precipitation(df_ppt, ppt_df_year, paths):
    
    # df_ppt = ppt_df
    
    df_ppt = df_ppt.reset_index()

    plt.figure(figsize = (10,5))
    
    p = sns.lineplot(data = df_ppt,
                     color='black',
                     x = 'month',
                     y = 'avg')
    
    plt.fill_between(x='month',y1='min',y2='max',data=df_ppt,facecolor='gray', alpha=0.2)
    
    plt.legend(labels=["Mean [2014-2023]","Max and Min [2014-2023]"], 
               borderaxespad=0, bbox_to_anchor=(0.99, 0.98), frameon=True, fontsize = 10)
    
    p.set_xlabel("Month", fontsize = 15)
    p.set_ylabel("Monthly Precipitation (mm)", fontsize = 15)
    plt.savefig(paths[0])
    
    
    plt.figure(figsize = (10,5))
    ax = sns.barplot(data=ppt_df_year, x="Ano", y="hourlyPrecipRate", color = 'c')
    
    ax.set_xlabel("Year", fontsize = 15)
    ax.set_ylabel("Precipitation (mm)", fontsize = 15)
    plt.savefig(paths[1])
    
    return

def multi_plot_temp(temp_df, temp_df_year, paths):
    
    temp_df.set_index('id', inplace=True)
    temp_df.sort_index(inplace=True)

    plt.figure(figsize = (10,5))
    
    p = sns.lineplot(data = temp_df,
                     color='black',
                     x = 'mes',
                     y = 'T_MEAN')
    
    plt.fill_between(x='mes',y1='T_MIN',y2='T_MAX',data=temp_df,facecolor='gray', alpha=0.2)
    
    plt.legend(labels=["Mean [2019-2023]","Max and Min [2019-2023]"], 
               borderaxespad=0, bbox_to_anchor=(0.99, 0.98), frameon=True, fontsize = 10)
    
    p.set_xlabel("Month", fontsize = 15)
    p.set_ylabel("Tempearture (C°)", fontsize = 15)
    plt.savefig(paths[0])
    
    
    plt.figure(figsize = (10,5))
    ax = sns.barplot(data=temp_df_year, x="ano", y="T_MEAN", color = 'c')
    
    ax.set_xlabel("Year", fontsize = 15)
    ax.set_ylabel("Temperature (C°)", fontsize = 15)
    plt.savefig(paths[1])
    
    return

multi_plot_precipitation(
    ppt_df, ppt_df_year,
        [r'C:\Projetos\bioflore\gbs-uganda\png\month_prec.png',
             r'C:\Projetos\bioflore\gbs-uganda\png\year_prec.png'])

multi_plot_temp(
    temp_df, temp_df_year,
                [r'C:\Projetos\bioflore\gbs-uganda\png\month_temp.png',
                     r'C:\Projetos\bioflore\gbs-uganda\png\year_temp.png'])
