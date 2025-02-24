# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 21:57:02 2024

@author: Gabriel
"""
import os
path = r"C:\Projetos\bioflore\bgci_II\spreadsheets"
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
from matplotlib.dates import DateFormatter 
#%%

ndvi = pd.read_csv("bspline_ndvi_addons.csv")
labels = list(set(ndvi.talhao))

r_list = []
g_list = []

# Loop over each label to process and plot the data
for label in labels:
    try:
        # Filter and pivot data
        dataframe = ndvi[ndvi.talhao == label]
        pivoted_data = dataframe.pivot(index='date', columns='talhao', values='ndvi_mean').reset_index()
        
        pivoted_data.date = pd.to_datetime(pivoted_data.date)
        pivoted_data.sort_values(by='date', inplace=True)
        
        date_non_na = pivoted_data[pivoted_data[label].notna()].iloc[0]['date']
        
        x = pivoted_data[label].ffill().bfill()
        y = pivoted_data['date']
        
        data = pd.DataFrame({'date': y, 'ndvi': x}, columns=['date', 'ndvi'])
        
        # Find troughs in the NDVI data
        troughs, _ = find_peaks(-data['ndvi'], distance=250, prominence=-0.000000000000000001)
    
        # Identify seasons based on trough indices
        seasons = [(data['date'].iloc[troughs[i-1]], data['date'].iloc[troughs[i]]) for i in range(1, len(troughs))]
        
        # Prepare for plotting
        fig, ax1 = plt.subplots(figsize=(14, 8))
        
        # Adjust font sizes
        plt.rcParams.update({'font.size': 14})  # Increase all fonts globally
        
        # Plot NDVI data
        ax1.set_xlabel('Date', fontsize=16)
        ax1.set_ylabel('NDVI', fontsize=16)
        ax1.plot(data['date'].iloc[troughs], data['ndvi'].iloc[troughs], 'o',
                 label='End of Season (EOS)', markersize=10, color='blue')
        ax1.plot(data['date'], data['ndvi'], label='NDVI', color='tab:blue', linestyle="--")
        ax1.tick_params(axis='y', labelsize=12)
        ax1.set_ylim([0, 1])
        ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
        plt.xticks(rotation=45, fontsize=12)
        plt.yticks(fontsize=12)
        
        # Secondary Y-axis for GPP
        ax2 = ax1.twinx()
        ax2.set_ylabel('Gross Primary Production (GPP)', fontsize=16)
        ax2.set_ylim([0, 500])
        
        integral_values = []
        season_centers = []
        
        for start, end in seasons:
            mask = (data['date'] >= start) & (data['date'] <= end)
            season_data = data.loc[mask]
            dates = season_data['date']
            ndvi_values = season_data['ndvi']
            
            # Integrate NDVI over the season
            integral_value = simps(ndvi_values.dropna(), dx=1)
            integral_values.append(integral_value)
            
            # Calculate the center date
            center_date = dates.iloc[len(dates) // 2]
            season_centers.append(center_date)
            
            # Shade the area under the curve
            ax1.fill_between(dates, ndvi_values, alpha=0.3)
        
        f = round(integral_values[0], 2)
        l = round(integral_values[-1], 2)
        r = round(((l / f) - 1) * 100, 2)
        
        r_list.append({'area': label, 'first_season': f, 'last_season': l, 'reason': r})
        g_list.append({'values': integral_values, 'area': label})
        
        # Plot integral values
        ax2.plot(season_centers, integral_values, 'o-', label='Season Integral',
                 color='tab:orange', linewidth=2.0, markersize=7)
        
        # Add title
        plt.title(f'NDVI Time Series with GPP values for {label.title()}', fontsize=18)
        
        # Add legends
        fig.legend(loc='upper left', bbox_to_anchor=(0.1, 0.9), fontsize=12)
        
        # Tight layout
        fig.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit title and legend
        
    except Exception as e:
        print(f"An error occurred for {label}: {e}")
        continue

# Save to CSV
df = pd.DataFrame(r_list)
df.to_csv("gpp_table.csv")  
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