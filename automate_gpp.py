# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 21:57:02 2024

@author: Gabriel
"""
import os
path = r"C:\Projetos\bioflore\ecosia\Ecosia"
os.chdir(path)

from datetime import datetime, timedelta
import pandas as pd
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

#%%
def plot_ndvi_curve(y, x):
    plt.plot([datetime.fromordinal(int(date)).date() for date in y], x, 'black', linewidth=3.5, label='Restoration Area')
    
    # Find and plot minimum points along the spline curve for each year within the range
    min_points_x = []
    min_points_y = []

    years = pd.to_datetime([datetime.fromordinal(int(date)).date() for date in y]).year
    unique_years = np.unique(years)
    
    for year in unique_years:
        in_year = years == year
        if np.sum(in_year) > 0:
            min_val_index = np.argmin(x[in_year])
            min_date = datetime.fromordinal(int(y[in_year][min_val_index]))
            min_val = x[in_year][min_val_index]
            
            min_points_x.append(min_date)
            min_points_y.append(min_val)
            
            # Plot the minimum point
            plt.plot(min_date, min_val, 'bo')
    
    # Connect the minimum points with lines
    if len(min_points_x) > 1:
        plt.plot(min_points_x, min_points_y, '#0101ff')
        
    return None


#%%
# x = np.linspace(0, 5, 1000)  # 100 pontos de 0 a 1
# # Gerar a função senoidal
# y = np.sin(2 * np.pi * x)  
# y = y+1

# data = DataFrame({'date': x, 'ndvi': y})

# # Detect troughs (which represent the lows in NDVI)
# troughs, _ = find_peaks(-data['ndvi'], distance=1, prominence=0.001)

# # Extracting start and end of seasons
# seasons = []
# for i in range(1, len(troughs)):
#     start_of_season = troughs[i-1]
#     end_of_season = troughs[i]
#     seasons.append((data['date'].iloc[start_of_season], data['date'].iloc[end_of_season]))

# # print("Identified Seasons (Start and End Dates):")
# # for season in seasons:
# #     print(f"Start: {season[0]}, End: {season[1]}")
    
# # Prepare for plotting
# fig, ax1 = plt.subplots(figsize=(10, 5))

# # Plot the NDVI data on the first Y-axis
# ax1.set_xlabel('Date')
# ax1.set_ylabel('NDVI')
# ax1.plot(data['date'].iloc[troughs], data['ndvi'].iloc[troughs], 'o',
#          label='End of Season (EOS)', markersize = 10, color = 'blue')
# ax1.plot(data['date'], data['ndvi'], label='NDVI', color='tab:blue',
#          linestyle = "--")
# ax1.tick_params(axis='y')
# # ax1.set_ylim([-1, 1])

# # Create a secondary Y-axis for the integral values
# # ax2 = ax1.twinx()
# # ax2.set_ylabel('Gross Primary Production (GPP)')
# # ax2.set_ylim([0, 500])
# integral_values = []
# season_centers = []

# for idx, (start, end) in enumerate(seasons, 1):
#     mask = (data['date'] >= start) & (data['date'] <= end)
#     season_data = data.loc[mask]  # Slicing the data for each season
#     dates = season_data['date']
#     ndvi_values = season_data['ndvi']
    
#     # Calculate the integral of the NDVI curve over this season
#     integral_value = simps(ndvi_values.dropna(), dx=1)  # Using the Simpson's rule for integration
    
#     integral_values.append(integral_value)
    
#     # Find the central date of this season
#     center_date = dates.iloc[len(dates) // 2]
#     season_centers.append(center_date)
    
#     # Optional: Shade the area under the NDVI curve for each season
#     ax1.fill_between(dates, ndvi_values, alpha=0.3)

# # Plot the integral values at the center dates on the secondary Y-axis
# # ax2.plot(season_centers, integral_values, 'o-', label='Season Integral',
# #          color='tab:orange', linewidth=2.0, markersize=7)
# # ax2.tick_params(axis='y')

# # Add legends
# fig.tight_layout()
# fig.legend(loc='upper left', bbox_to_anchor=(0.1,0.9))

# plt.title('Senoidal test')
# # plt.savefig(r"C:\Projetos\bioflore\ecosia\plots_gpp"+ f"//gpp_{label}.png")

#%%

ndvi = pd.read_csv(path+r"/bspline_ndvi.csv")
labels = list(set(ndvi.talhao))

# labels = [i for i in labels if i.startswith('ESP_')]
r_list = []

for label in labels:
    try:
        # label = "ESP_05"
        dataframe = ndvi[ndvi.talhao==label]
        # data.rename(columns={'ndvi_mean':'ndvi'}, inplace=True)
        pivoted_data = dataframe.pivot(index='date', columns='talhao', values='ndvi_mean').reset_index()
        
        pivoted_data.date = pd.to_datetime(pivoted_data.date, format='%Y-%m-%d')
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
        troughs, _ = find_peaks(-data['ndvi'], distance=250, prominence=0.01)
    
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
        # Plot the integral values at the center dates on the secondary Y-axis
        ax2.plot(season_centers, integral_values, 'o-', label='Season Integral',
                  color='tab:orange', linewidth=2.0, markersize=7)
        ax2.tick_params(axis='y')
        
        # Add legends
        fig.tight_layout()
        fig.legend(loc='upper left', bbox_to_anchor=(0.1,0.9))
        
        plt.title(f'NDVI Time Series with GPP value {label}')
        plt.savefig(r"C:\Projetos\bioflore\ecosia\plots_gpp_t"+ f"//gpp_{label}.png")
        plt.show()
    except:
        continue
   
df = pd.DataFrame(r_list)
df.to_csv(r"C:\Projetos\bioflore\ecosia\plots_gpp_t\gpp_table.csv")  
