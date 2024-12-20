# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:18:49 2024

@author: Gabriel
"""
import geopandas as gpd
import io
import os
import zipfile
import ee
import requests
from shapely.geometry import mapping
from google.oauth2 import service_account
from shapely.geometry import Polygon, MultiPolygon
from shapely.geometry import box
import numpy as np

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

def search_on_gee(extent, source, band, name,
                  mosaic=False, buffer_km=0, scale =10):
    
    # extent=extent
    # source="WorldPop/GP/100m/pop_age_sex_cons_unadj"
    # band='population'
    # name= name
    # mosaic=True
    
    if buffer_km > 0:
        buffer_distance = buffer_km/ 111
        extent['geometry'] = extent.geometry.buffer(buffer_distance)
    
    aoi_2d = remove_z(extent)
    polygon = geometry_to_ee(aoi_2d)
    
    if mosaic:
        dataset = (ee.ImageCollection(source)
                   .select(band)     
                 .filterBounds(polygon).mosaic())
    else:
        dataset = (ee.Image(source)
                   .select(band))

    file_name = f'{band}.tif'
    directory = f'/tmp/{band}/{name}'

    if not os.path.exists(directory):
        os.makedirs(directory)

    output_dir = f'{directory}/{file_name}'
    
    image = dataset.select(band)

    print(image.getInfo())

    img_crs =  image.projection().crs().getInfo()
    
    url = image.getDownloadURL({
        'scale': scale,
        'crs': img_crs,
        'region': polygon.toGeoJSONString(),
    })

    response = requests.get(url, stream=True, timeout=1000000)

    with zipfile.ZipFile(io.BytesIO(response.content)) as collection_zip:
        collection_zip.extractall(output_dir)
        # shutil.rmtree(f"{directory}/extracted")
    return output_dir

def get_extent_as_geodf(gdf):
    # Get the total bounds of the GeoDataFrame
    bounds = gdf.total_bounds  # Returns an array: [minx, miny, maxx, maxy]

    # Create a bounding box (rectangle) using shapely
    bbox = box(bounds[0], bounds[1], bounds[2], bounds[3])

    # Create a new GeoDataFrame containing just this bounding box
    extent_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs=gdf.crs)

    return extent_gdf

def get_earth_engine_connection(key: str):
    credentials = service_account.Credentials.from_service_account_file(key)
    scoped_credentials = credentials.with_scopes(
        ['https://www.googleapis.com/auth/cloud-platform'])
    return scoped_credentials

def split_aoi(aoi_gdf, grid_size=0.5):
    # Convert AOI to WGS84 (EPSG:4326)
    aoi_gdf = aoi_gdf.to_crs('EPSG:4326')

    # Get bounds of the AOI
    minx, miny, maxx, maxy = aoi_gdf.total_bounds

    # Create grids
    x_coords = np.arange(minx, maxx, grid_size)
    y_coords = np.arange(miny, maxy, grid_size)
    grid_polygons = [
        box(x, y, x + grid_size, y + grid_size)
        for x in x_coords
        for y in y_coords
    ]

    grid_gdf = gpd.GeoDataFrame({'geometry': grid_polygons}, crs='EPSG:4326')

    # Intersect the grid with AOI
    intersected_gdf = gpd.overlay(grid_gdf, aoi_gdf, how='intersection')
    return intersected_gdf

gdf = gpd.read_file(r"C:\Projetos\bioflore\bgci_II\geo\sites\all_sites\all_sites.shp")
name = 'bgciII'

KEY = r"C:\Projetos\bioflore\restoration_monitoring\bioflore-ee.json"
session = get_earth_engine_connection(KEY)
ee.Initialize(session)

extent = get_extent_as_geodf(gdf.to_crs(4326))

dem = search_on_gee(extent, 'NASA/NASADEM_HGT/001', 'elevation',
                    name, mosaic=False, scale = 30)
soil_020 = search_on_gee(extent, "ISDASOIL/Africa/v1/texture_class",
                         'texture_0_20', name, mosaic=False, scale=30)
soil_2050 = search_on_gee(extent, "ISDASOIL/Africa/v1/texture_class",
                         'texture_20_50', name, mosaic=False, scale=30)

splited_aoi = split_aoi(extent, 0.2)

for idx, row in splited_aoi.iterrows():
    aoi = gpd.GeoDataFrame({'geometry': [row.geometry]}, crs=splited_aoi.crs)
    try:
        output_dir = search_on_gee(
            extent=aoi,
            source= "ESA/WorldCover/v200",
            band= 'Map',
            name=f"{name}_{idx}",  # Generate unique name for each operation
            mosaic=True,
            scale= 10
        )
        print(f"Processed extent {idx}, output saved to {output_dir}")
    except Exception as e:
        print(f"Error processing extent {idx}: {e}")
   
data_tuples = [("ISDASOIL/Africa/v1/texture_class",'texture_0_20'),
             ("ISDASOIL/Africa/v1/texture_class", 'texture_20_50'),
             ("ISDASOIL/Africa/v1/bedrock_depth", 'mean_0_200')]

for idx, row in gdf.iterrows():
    aoi = get_extent_as_geodf(
        gpd.GeoDataFrame({'geometry': [row.geometry]}, crs=gdf.crs))
    label = row.Name
    for tup in data_tuples:
        try:
            output_dir = search_on_gee(
                extent=aoi,
                source= tup[0],
                band=   tup[1],
                name=f"{name}_{label}",  # Generate unique name for each operation
                mosaic=False,
                scale= 30
            )
            print(f"Processed extent {idx}, output saved to {output_dir}")
        except Exception as e:
            print(f"Error processing extent {idx}: {e}")

