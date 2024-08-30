# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 11:58:13 2024

@author: Gabriel
"""
import pandas as pd
import rasterio
from rasterio import mask
import numpy as np
from datetime import datetime
import pyproj
from shapely.ops import transform

def get_ndvi_mean(
        fname: str,
        shape: pd.core.series.Series,
        aoi_crs: pyproj.crs.crs.CRS)-> tuple:
    
    
    # path = 'C:\\Projetos\\bioflore\\ecosia\\PSdata\\'
    src1 = rasterio.open(fname)
    project = pyproj.Transformer.from_crs(aoi_crs, src1.crs, always_xy=True).transform
    reprojected = transform(project, shape['geometry'])
    
    
    geom = [reprojected.__geo_interface__]
    clipped_array, clipped_transform = mask(dataset=src1, shapes=geom,
                                            crop=True)

    sr_data = clipped_array.read()

    src2 = rasterio.open(fname.replace('.tif', '_QA.tif'))
    qa_data = src2.read(1)
    
    mask_array = np.where(qa_data == 0, 0, 1)
    
    # Apply the mask to the surface reflectance raster
    masked_sr_data = sr_data * mask_array
    
    nir_band = masked_sr_data[3]
    red_band = masked_sr_data[2]
    ndvi = (nir_band - red_band) / (nir_band + red_band)
    ndvi_mean = np.nanmean(ndvi, dtype=np.float64)
    
    dn = fname.split('_')[-2] + fname.split('_')[-1].replace('.tif','')
    date = datetime.strptime(dn, '%Y%m%d%H%M%S')
    talhao = fname[:6]
    
    result = (date, ndvi_mean, talhao)
    
    # os.remove(directory +'/'+ img_name)
    # os.remove(directory +'/'+ qa_name)
    
    return result

