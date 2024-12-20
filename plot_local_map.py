# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:07:50 2024

@author: Gabriel
"""
# import os
# import geopandas as gpd
# import contextily as ctx
# import matplotlib.pyplot as plt
# from shapely.geometry import box

# os.chdir(r"C:\Projetos\bioflore\gbs-kenya")
# aoi = gpd.read_file("geo/shp/restoration_sites_ms.shp")
# aoi = aoi.to_crs(epsg=3857)

# for idx, row in aoi.iterrows():
#     polygons = gpd.GeoDataFrame([row], crs=aoi.crs)
    
#     xmin, ymin, xmax, ymax = polygons.total_bounds
#     bbox = box(xmin, ymin, xmax, ymax)
#     bbox_gdf = gpd.GeoSeries([bbox], crs=polygons.crs)
#     buffer_distance = 1000  # Ensure this is 5 km
#     buffered_bbox = bbox_gdf.buffer(buffer_distance)

#     fig, ax = plt.subplots(figsize=(10, 10))

#     polygons['label'] = polygons['label'].str.title()
#     polygons.plot(ax=ax, edgecolor='black', alpha=1, column='label', legend=True)
#     buffered_bbox.boundary.plot(ax=ax, edgecolor='none', alpha=0)

#     # Add label within polygon
#     centroid = row['geometry'].centroid
#     ax.text(centroid.x, centroid.y, row['label'], fontsize=10,
#             ha='center', va='center', color='black', weight='bold')

#     # Add basemap
#     ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

#     buffer_total_bounds = buffered_bbox.total_bounds
#     ax.set_xlim(buffer_total_bounds[0], buffer_total_bounds[2])
#     ax.set_ylim(buffer_total_bounds[1], buffer_total_bounds[3])

#     # Hide x and y ticks
#     ax.set_xticks([])
#     ax.set_yticks([])

#     # Add a custom scale bar
#     scalebar_text = '1 km'  # Label for the scale bar
#     scalebar_length = 1000  # Length in meters

#     # Coordinates for the scalebar, you might need to adjust these
#     scalebar_x = buffer_total_bounds[0] + 0.01 * (buffer_total_bounds[2] - buffer_total_bounds[0])
#     scalebar_y = buffer_total_bounds[1] + 0.05 * (buffer_total_bounds[3] - buffer_total_bounds[1])

#     ax.plot([scalebar_x, scalebar_x + scalebar_length], [scalebar_y, scalebar_y], color='black', lw=3)
#     ax.text(scalebar_x + scalebar_length / 2, scalebar_y + scalebar_length * 0.02, scalebar_text,
#             horizontalalignment='center', verticalalignment='bottom')

#     plt.legend(bbox_to_anchor=(0.5, 1.05),
#                fancybox=True, shadow=True)

#     plt.show()
#     fig.savefig(f'png/local_{idx}.png')
    
    
import os
import geopandas as gpd
import contextily as ctx
import matplotlib.pyplot as plt
from shapely.geometry import box

os.chdir(r"C:\Projetos\bioflore\gbs-kenya")
aoi = gpd.read_file("geo/shp/restoration_sites.shp")
aoi = aoi.to_crs(epsg=3857)

# Determine the number of polygons and define grid size
num_polygons = len(aoi)
num_cols = 5  # Number of columns in the grid
num_rows = (num_polygons + num_cols - 1) // num_cols  # Calculate rows needed

fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 3))
axes = axes.flatten()  # Flatten axes array for easy iteration

for ax, (idx, row) in zip(axes, aoi.iterrows()):
    polygons = gpd.GeoDataFrame([row], crs=aoi.crs)
    
    # Calculate bounds for the current polygon
    xmin, ymin, xmax, ymax = polygons.total_bounds
    bbox = box(xmin, ymin, xmax, ymax)
    bbox_gdf = gpd.GeoSeries([bbox], crs=polygons.crs)
    buffer_distance = 1000  # Ensure this is 5 km
    buffered_bbox = bbox_gdf.buffer(buffer_distance)

    # Plotting each polygon
    polygons['label'] = polygons['label'].str.title()
    polygons.plot(ax=ax, edgecolor='black', alpha=1, column='label', legend=False)
    buffered_bbox.boundary.plot(ax=ax, edgecolor='none', alpha=0)

    # Add label within the polygon
    centroid = row['geometry'].centroid
    ax.text(centroid.x, (centroid.y-500), row['label'], fontsize=10,
            ha='center', va='center', color='black', weight='bold')

    # Add basemap
    ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

    buffer_total_bounds = buffered_bbox.total_bounds
    ax.set_xlim(buffer_total_bounds[0], buffer_total_bounds[2])
    ax.set_ylim(buffer_total_bounds[1], buffer_total_bounds[3])

    # Hide x and y ticks for clarity
    ax.set_xticks([])
    ax.set_yticks([])

# Hide any empty subplots
for j in range(idx + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.show()
fig.savefig('png/combined_polygons.png')
