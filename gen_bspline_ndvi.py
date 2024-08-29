# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:25:16 2024

@author: Gabriel
"""

import pandas as pd
import numpy as np
from scipy.interpolate import BSpline
from patsy import dmatrix
import statsmodels.api as sm

pivoted_data  = pd.read_csv(r"C:\Projetos\bioflore\ecosia\i.txt")

label = 'ndvi_mean'

# Ensure dates are in datetime format
pivoted_data['date'] = pd.to_datetime(pivoted_data['date'])

# Sort the DataFrame by date
pivoted_data = pivoted_data.sort_values(by='date')

# Find the first non-NA date for the given label
date_non_na = pivoted_data.loc[~pivoted_data[label].isna(), 'date'].iloc[0]

# Define the specific range for the x-axis
start_date = pivoted_data.iloc[0]['date']
end_date = pivoted_data.iloc[-1]['date']

# Forward fill interpolation
pivoted_data[label] = pivoted_data[label].fillna(method='ffill')

# Prepare data for spline fitting
x = pivoted_data[label]
y = pivoted_data['date']

# Convert date to numeric for spline fitting
y_numeric = (y - y.min()).dt.days

# Degree and degrees of freedom
degree = 4
df = 15

# Fit B-spline
spline_basis = dmatrix("bs(train, degree=degree, df=df, include_intercept=False)", {"train": y_numeric}, return_type='dataframe')
fit_bs = sm.OLS(x, spline_basis).fit()

# Predictions for a smooth curve
y_pred = pd.date_range(start=date_non_na, end=end_date, freq='D')
y_pred_numeric = (y_pred - y.min()).days
spline_basis_pred = dmatrix("bs(train, degree=degree, df=df, include_intercept=False)", {"train": y_pred_numeric}, return_type='dataframe')
x_pred = fit_bs.predict(spline_basis_pred)

# Creating the resulting DataFrame
l = [label] * len(x_pred)
result_df = pd.DataFrame({'date': y_pred, 'ndvi_mean': x_pred, 'talhao': l})

import matplotlib.pyplot as plt
import seaborn as sns

# Ensure seaborn style for better aesthetics
sns.set()

# Plotting the result_df
plt.figure(figsize=(12, 6))

# Plot the original data
plt.plot(pivoted_data['date'], pivoted_data[label], 'o', label='Original data', markersize=5)

# Plot the smoothed curve
plt.plot(result_df['date'], result_df['ndvi_mean'], label='Smoothed B-spline fit', color='red', linewidth=2)

# Adding labels, title, legend, etc.
plt.xlabel('Date')
plt.ylabel('NDVI Mean')
plt.title('NDVI Mean Over Time with B-Spline Fit')
plt.legend()
plt.grid(True)

# Show plot
plt.show()
