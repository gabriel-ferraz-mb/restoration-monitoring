# Here's a step-by-step guide to set `R_HOME` and make sure `rpy2` works correctly:

# ### Set `R_HOME` (Windows):

# 1. **Find R Installation Path**:
#     Locate the directory where R is installed. The path typically looks something like `C:\Program Files\R\R-4.1.0`.

# 2. **Set `R_HOME` Environment Variable**:
#     - **Manually via Windows Settings**:
#         1. Right-click on "This PC" or "My Computer" on the Desktop or in File Explorer.
#         2. Select "Properties".
#         3. Click on "Advanced system settings".
#         4. Click on the "Environment Variables" button.
#         5. In the "System variables" section, click "New".
#         6. Enter `R_HOME` as the variable name and your R installation path (e.g., `C:\Program Files\R\R-4.1.0`) as the variable value.
#         7. Click OK to close all windows.

#     - **Programmatically in Python**:
#         You can set `R_HOME` within your Python script before importing `rpy2`:
#         ```python
#         import os

#         # Make sure to replace the path with your actual R installation path
#         os.environ['R_HOME'] = 'C:\\Program Files\\R\\R-4.1.0'

#         import rpy2.robjects as robjects
#         from rpy2.robjects import pandas2ri
#         from rpy2.robjects.packages import importr
#         from rpy2.robjects.vectors import StrVector
        
#         pandas2ri.activate()

#         # Rest of your code follows...
#         ```

# ### Complete Python Script with `R_HOME` Set and `rpy2` Integration:

# ```python
# Set R_HOME environment variable
import os

# Make sure to replace the path with your actual R installation path
os.environ['R_HOME'] = "C:\Program Files\R\R-4.3.3"

# Import necessary libraries
# Import necessary libraries
import rasterio
import numpy as np
import pandas as pd
from os import listdir
from os.path import join
from rasterio.plot import show
import matplotlib.pyplot as plt
from datetime import datetime
from statsmodels.nonparametric.smoothers_lowess import lowess
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

# Activate the conversion from pandas to R dataframes
pandas2ri.activate()

# Import necessary R packages
r_base = importr('base')
raster = importr('raster')
RStoolbox = importr('RStoolbox')
terra = importr('terra')

def write_raster(data, file_path, ref_file):
    with rasterio.open(ref_file) as src:
        profile = src.profile
        profile.update(dtype=rasterio.float32, count=1)
        with rasterio.open(file_path, 'w', **profile) as dst:
            dst.write(data, 1)
            
def calc_ndfi(name, em):
    product = name.replace("PSdata", "NDFI").replace(".tif", "_NDFI_test.tif")
    
    if os.path.exists(product):
        return None
    
    with rasterio.open(name) as src:
        data = src.read()
        height, width = src.shape
        # show(src.read([3,2,1]))
    
    m_file = name.replace(".tif", "_QA.tif")
    
    with rasterio.open(m_file) as src:
        m = src.read()
  
    
    imgBrick = data * m
    contains_non_zero = np.any(imgBrick != 0, axis=(1, 2))
    if not contains_non_zero:
        return None
    
     # Convert the data to R objects
    r_name = robjects.StrVector([name])
    r_data = raster.brick(r_name)
    r_endmembers = pandas2ri.py2rpy(em)
    
    # Apply the mesma function in R
    imgFracao = RStoolbox.mesma(r_data, r_endmembers)
    imgFracao_layers = terra.as_list(imgFracao)
    
    veg = np.array(terra.values(imgFracao_layers[0]).reshape(height, width))
    agua = np.array(terra.values(imgFracao_layers[1]).reshape(height, width))
    palha = np.array(terra.values(imgFracao_layers[2]).reshape(height, width))
    areia = np.array(terra.values(imgFracao_layers[3]).reshape(height, width))
    
    gv_shade = veg / (1 - agua)
    NDFI = (gv_shade - (palha + areia)) / (palha + areia + gv_shade)
    write_raster(NDFI, product, name)
    return None
   # Write NDFI as a raster

# Set working directory and list files
work_dir = 'C:\\Projetos\\bioflore\\ecosia'
os.chdir(work_dir)

files = ['C:\\Projetos\\bioflore\\ecosia\\PSdata\\' + f for f in listdir("PSdata") if 'ESP' in f and 'QA.tif' not in f]

# Read and process endmembers
endmembers = pd.read_csv("em.txt")
endmembers.columns = ["class", "Blue", "Green", "Red", "NIR"]
endmembers.set_index('class', inplace = True)

for f in files:
    calc_ndfi(f, endmembers)
    


    # Add additional plotting logic here if needed

# Ensure you handle specific operations not directly covered by basic libraries.
# ```

# By setting the `R_HOME` environment variable in the script, you ensure that `rpy2` can locate the R installation without requiring manual configuration each time the script is run. This setup should work properly assuming R is installed and accessible at the specified path.