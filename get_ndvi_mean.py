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
import os
from os import listdir

def get_ndvi_mean(fname):
    
    path = 'C:\\Projetos\\bioflore\\ecosia\\PSdata\\'
    src1 = rasterio.open(path + fname)
    sr_data = src1.read()

    src2 = rasterio.open(path + fname.replace('.tif', '_QA.tif'))
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

work_dir = 'C:\\Projetos\\bioflore\\ecosia'
os.chdir(work_dir)

files = [ f for f in listdir("PSdata")  if 'QA.tif' not in f]

result_list = []
for f in files:
    result_list.append(get_ndvi_mean(f))

df = pd.DataFrame(result_list, columns=['date', 'ndvi_mean', 'talhao'])

df.to_csv(work_dir + '/' + 'resultado_22062024.csv', index=False)