# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:07:50 2024

@author: Gabriel
"""
import os
import geopandas as gpd
import contextily as ctx
import matplotlib.pyplot as plt
from shapely.geometry import box

os.chdir(r"C:\Projetos\bioflore\gbs-uganda")
polygons = gpd.read_file("geo/shp/restoration_sites.shp")
polygons = polygons.to_crs(epsg=3857)

xmin, ymin, xmax, ymax = polygons.total_bounds
bbox = box(xmin, ymin, xmax, ymax)
bbox_gdf = gpd.GeoSeries([bbox], crs=polygons.crs)
buffer_distance = 5000  # Ensure this is 5 km
buffered_bbox = bbox_gdf.buffer(buffer_distance)

fig, ax = plt.subplots(figsize=(10, 10))

polygons['Name'] = polygons['Name'].str.title()
polygons.plot(ax=ax, edgecolor='red', alpha=0.5, column='Name', legend=True)
buffered_bbox.boundary.plot(ax=ax, edgecolor='none', alpha=0)

ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

buffer_total_bounds = buffered_bbox.total_bounds
ax.set_xlim(buffer_total_bounds[0], buffer_total_bounds[2])
ax.set_ylim(buffer_total_bounds[1], buffer_total_bounds[3])

# Hide x and y ticks
ax.set_xticks([])
ax.set_yticks([])

# Add a custom scale bar
scalebar_text = '1 km'  # Label for the scale bar
scalebar_length = 1000  # Length in meters

# Coordinates for the scalebar, you might need to adjust these
scalebar_x = buffer_total_bounds[0] + 0.01 * (buffer_total_bounds[2] - buffer_total_bounds[0])
scalebar_y = buffer_total_bounds[1] + 0.05 * (buffer_total_bounds[3] - buffer_total_bounds[1])

ax.plot([scalebar_x, scalebar_x + scalebar_length], [scalebar_y, scalebar_y], color='black', lw=3)
ax.text(scalebar_x + scalebar_length/2, scalebar_y + scalebar_length * 0.02, scalebar_text, horizontalalignment='center', verticalalignment='bottom')

plt.show()
fig.savefig('png/local.png')
