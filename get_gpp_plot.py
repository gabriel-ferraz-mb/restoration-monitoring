# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 21:57:02 2024

@author: Gabriel
"""
import os
path = r"C:\Projetos\bioflore\gbs-kenya\spreadsheets"
os.chdir(path)

from datetime import datetime, timedelta
import pandas as pd
import geopandas as gpd
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
from numpy import trapz
from scipy.interpolate import BSpline, make_interp_spline
from pandas.tseries.offsets import DateOffset
from os import listdir
from scipy.signal import find_peaks
from scipy.integrate import simps
from dateutil.relativedelta import relativedelta
import seaborn as sns
#%%

ndvi = pd.read_csv("bspline_ndvi_lmite.csv")
labels = list(set(ndvi.talhao))

shp= gpd.read_file("../geo/shp/restoration_sites.shp")
r_list = []

g_list = []

for label in labels:
    try:
        # label = labels[0]
        dataframe = ndvi[ndvi.talhao==label]
        
        pivoted_data = dataframe.pivot(index='date', columns='talhao', values='ndvi_mean').reset_index()
        
        pivoted_data.date = pd.to_datetime(pivoted_data.date)
        pivoted_data.sort_values(by='date', inplace = True)
        
        date_non_na = pivoted_data[pivoted_data[label].notna()].iloc[0]['date']
        start_date = pivoted_data['date'].min()
        end_date = pivoted_data['date'].max()
        
        x = pivoted_data[label].ffill().bfill()
        y = pivoted_data['date']
        
        data = pd.DataFrame({'date': y, 'ndvi': x}, columns=['date', 'ndvi'])
        
        # Detect peaks (which represent the highs in NDVI)
        # peaks, _ = find_peaks(data['ndvi'], distance=300, prominence=[-0.1,1])
        
        # Detect troughs (which represent the lows in NDVI)
        troughs, _ = find_peaks(-data['ndvi'], distance=250, prominence=-0.000000000000000001)
    
        # Extracting start and end of seasons
        seasons = []
        for i in range(1, len(troughs)):
            start_of_season = troughs[i-1]
            end_of_season = troughs[i]
            seasons.append((data['date'].iloc[start_of_season], data['date'].iloc[end_of_season]))
        
        print("Identified Seasons (Start and End Dates):")
        for season in seasons:
            print(f"Start: {season[0]}, End: {season[1]}")
            
        # Prepare for plotting
        fig, ax1 = plt.subplots(figsize=(14, 8))
        
        # Plot the NDVI data on the first Y-axis
        ax1.set_xlabel('Date')
        ax1.set_ylabel('NDVI')
        ax1.plot(data['date'].iloc[troughs], data['ndvi'].iloc[troughs], 'o',
                  label='End of Season (EOS)', markersize = 10, color = 'blue')
        ax1.plot(data['date'], data['ndvi'], label='NDVI', color='tab:blue',
                  linestyle = "--")
        ax1.tick_params(axis='y')
        ax1.set_ylim([0, 1])
        # if 'reference' not in label:
        #     planting_date = datetime.strptime(shp[shp['identifier_akvo'] == label]['planting_date'].values[0],
        #                                 '%Y-%m-%d')
        #     ax1.axvline(x=planting_date, color='red', linestyle='--', label='Planting Date') 
        
        # Create a secondary Y-axis for the integral values
        ax2 = ax1.twinx()
        ax2.set_ylabel('Gross Primary Production (GPP)')
        ax2.set_ylim([0, 500])
        integral_values = []
        season_centers = []
        
        for idx, (start, end) in enumerate(seasons, 1):
            mask = (data['date'] >= start) & (data['date'] <= end)
            season_data = data.loc[mask]  # Slicing the data for each season
            dates = season_data['date']
            ndvi_values = season_data['ndvi']
            
            # Calculate the integral of the NDVI curve over this season
            integral_value = simps(ndvi_values.dropna(), dx=1)  # Using the Simpson's rule for integration
            
            integral_values.append(integral_value)
            
            # Find the central date of this season
            center_date = dates.iloc[len(dates) // 2]
            season_centers.append(center_date)
            
            # Optional: Shade the area under the NDVI curve for each season
            ax1.fill_between(dates, ndvi_values, alpha=0.3)
        f =  round(integral_values[0],2)
        l = round(integral_values[-1],2)
        r = round(((l/f)-1)*100,2)
        p = {'area':label, 'first_season':f,
             'last_season':l , 'reason':r}
        r_list.append(p)
        g_list.append({'values':integral_values, 'area':label})
        # Plot the integral values at the center dates on the secondary Y-axis
        ax2.plot(season_centers, integral_values, 'o-', label='Season Integral',
                  color='tab:orange', linewidth=2.0, markersize=7)
        # ax2.annotate("GPP value: " + str(round(integral_values[0],2)),
        #              xy=(start + relativedelta(months=2),integral_values[0]-50),
        #              weight='bold',
        #              size = 14)
        # ax2.tick_params(axis='y')
        
        # Add legends
        fig.tight_layout()
        fig.legend(loc='lower right', bbox_to_anchor=(0.9,0.1))
        
        plt.title(f'NDVI Time Series with GPP value {label.title()}')
        # plt.show()
        # plt.savefig(f"gpp_{label}.png")
    except:
        continue
   
    
df = pd.DataFrame(r_list)
df.to_csv(path + "//gpp_table.csv")  
#%%
data_expanded = []
for entry in g_list:
    for i, value in enumerate(entry['values']):
        data_expanded.append({'Time': f'Season {i + 1}', 'GPP': value, 'Area': entry['area']})

df = pd.DataFrame(data_expanded)
df['Area'] = df['Area'].str.title()
# Plotting
plt.figure(figsize=(10, 6))
sns.set(style="darkgrid")

sns.lineplot(x='Time', y='GPP', hue='Area', 
             marker='o', data=df)

plt.title('GPP along time series by area')
plt.xlabel('Time')
plt.ylabel('GPP')
plt.xticks(rotation=45)
plt.legend(title='Area', bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust the layout
plt.tight_layout()
plt.show()