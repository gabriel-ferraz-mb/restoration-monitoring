import io
import os
import shutil
import zipfile
import numpy as np
import ee
import pandas as pd
import geopandas as gpd
import rasterio
import requests
from datetime import datetime
from shapely.geometry import mapping
from google.oauth2 import service_account
from shapely.geometry import Polygon, MultiPolygon

def remove_z(geometry):
    if geometry.is_empty:
        return geometry  # Handle empty geometries if necessary

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

def get_earth_engine_connection(key: str):
    credentials = service_account.Credentials.from_service_account_file(key)
    scoped_credentials = credentials.with_scopes(
        ['https://www.googleapis.com/auth/cloud-platform'])
    return scoped_credentials

def get_bands(band_list: list):
    bands = []
    bands.append({
            "code": "B2",
            "band": "Blue",
            "resolution": "10m",
    })
    bands.append({
            "code": "B3",
            "band": "Green",
            "resolution": "10m",
    })
    bands.append({
            "code": "B4",
            "band": "Red",
            "resolution": "10m",
    })
    bands.append({
            "code": "B8",
            "band": "NIR",
            "resolution": "10m",
    })
    bands.append({
            "code": "B11",
            "band": "SWIR16",
            "resolution": "20m",
    })
    bands.append({
            "code": "B12",
            "band": "SWIR22",
            "resolution": "20m",
    })
    bands.append({
            "code": "SCL",
            "band": "Scene Classification Map",
            "resolution": "20m",
    })

    if len(band_list) == 0:
        return bands
    return [band for band in bands if band['code'] in band_list]

def sentinel2_get_images(
    aoi: gpd.GeoDataFrame,
    start_date: str,
    end_date: str, 
    granularity: str,
    name: str,
    band_list: list = [],
    cloud_cover: int = 101,
    KEY: str = r"C:\Projetos\bioflore\gbs-madagascar\bioflore-ee.json",
    collection: str = "COPERNICUS/S2_SR_HARMONIZED"
):
    
    session = get_earth_engine_connection(KEY)
    ee.Initialize(session)
    
    # aoi_2d = remove_3th_dimension(aoi.to_crs(4326))
    
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

    dataset = (ee.ImageCollection(collection)
            .filterDate(start_date, end_date)
            .filterBounds(polygon)
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_cover)))
    

    if len(band_list) > 0:
        dataset = dataset.select(band_list)

    if granularity == "month":
        dates = pd.date_range(start=start_date, end=end_date, freq="MS")
    elif granularity == "year":
        dates = pd.date_range(start=start_date, end=end_date, freq="YS")

    image_list = ee.List([])
    for date in dates:
        
        year = date.year
        month = date.month
        
        if granularity == "month":
            start_range = f"{year}-{month:02d}-01"
            end_range = (date + pd.offsets.MonthEnd()).strftime("%Y-%m-%d")
        elif granularity == "year":
            start_range = f"{year}-01-01"
            end_range = f"{year}-12-31"
            
        filtered_dataset = dataset.filterDate(start_range, end_range)
        sorted_by_cloud_cover = filtered_dataset.sort('CLOUDY_PIXEL_PERCENTAGE')
        least_cloudy_image = sorted_by_cloud_cover.first()
        try:
            if least_cloudy_image.getInfo()['type'] == 'Image':
                image_list = image_list.add(least_cloudy_image)
        except:
            continue

    if image_list.size().getInfo() == 0:
        image_list = dataset.toList(1, 0)

    description = []
    # i = 0
    # while image_list.size().getInfo() > 0:
    for i in range(0,image_list.size().getInfo()):
        image = ee.Image(image_list.get(i)).clip(polygon)
        date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo()
        file_name = f'{name}_{str(date)}.tif'
        directory = f'/tmp/sentinel2/{name}'

        if not os.path.exists(directory):
            os.makedirs(directory)

        output_file = f'{directory}/{file_name}'
        
        if len(band_list) > 0:
            img_crs = image.select(band_list[0]).projection().crs().getInfo()
        else:
            img_crs = image.select('B2').projection().crs().getInfo()

        url = image.getDownloadURL({
            'scale': 10,
            'crs': img_crs,
            'region': polygon.toGeoJSONString(),
        })

        response = requests.get(url, stream=True, timeout=1000000)

        with zipfile.ZipFile(io.BytesIO(response.content)) as collection_zip:
            collection_zip.extractall(f"{directory}/extracted")
            extracted_files = collection_zip.namelist()
            # merge all bands into a single raster using rasterio
            band_files = [f'{directory}/extracted/{band}' for band in extracted_files]
            output_file = f'{directory}/{file_name}'
            
            with rasterio.open(band_files[0]) as seed:
                profile = seed.profile
    
                # Atualizar o profile para multi-banda
                profile.update(count=len(band_files))
                
                # Criar uma matriz 3D para armazenar todas as bandas
                all_bands = np.zeros((len(band_files), seed.height, seed.width), dtype=seed.dtypes[0])
                
                # Ler cada banda e armazenar na matriz 3D
                for idx, file in enumerate(band_files):
                    with rasterio.open(file) as src:
                        all_bands[idx, :, :] = src.read(1)
            
            with rasterio.open(output_file, 'w', **profile) as dst:
                dst.write(all_bands)

            shutil.rmtree(f"{directory}/extracted")

        description.append(
            {
                'name': file_name,
                'path': output_file,
                'date': [date],
                'granularity': granularity,
                'satellite': 'sentinel 2',
                'bands': get_bands(band_list)
            }
        )
        # i += 1
        # image_list = dataset.toList(1, i)

    return description

path = r"C:\Projetos\bioflore\gbs-uganda\geo\shp\restoration_sites.shp"
aoi = gpd.read_file(path)
aoi['geometry'] = aoi['geometry'].apply(remove_z)
aoi.to_crs(4326, inplace=True)

# aoi =aoi[aoi['Name']== 'degraded'].reset_index()

for i in range(0,len(aoi)):
    site =aoi.iloc[[i]]
    name = site['Name'][i]
    # sd = site['Data'][i]
    sd = '2018-01-01'
    ed = datetime.today().strftime('%Y-%m-%d')
    
    if os.path.isdir(f'C:\\tmp\\sentinel2\\{name}'):
        print (f'{name} already exists')
        continue
    else:
        print(f"Starting download site {name}")
    
    try:
        d = sentinel2_get_images(
            aoi = site,
            start_date = sd,
            end_date = ed, 
            granularity = 'month',
            name= name,
            band_list = ['B2','B3','B4','B8','B11','B12','SCL'],
            cloud_cover = 100,
            KEY = r"C:\Projetos\bioflore\restoration_monitoring\bioflore-ee.json",
            collection = "COPERNICUS/S2_SR_HARMONIZED"
        )
    except Exception as e:
        print(f"Site {name} failed due to: {e}")
        continue






