"""
California Industrial and Auto Pollution API
Author: Florin Langer, UC Berkeley Renewable and Appropriate Energy Lab with Lawrence Berkeley National Laboratory
Sources: California Air Resources Board
"""

# %matplotlib inline # For Jupyter Notebook
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from scipy import ndimage
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
pylab.rcParams['figure.figsize'] = 8, 6
import folium
from folium.plugins import MarkerCluster
from folium import IFrame
import shapely
import unicodedata
import pysal as ps
import json
import os.path
import math
import ogr
import shapefile as shp # pyshp
import math
from helpers import save_name, county_code, read_census_shape, read_facilities,\
agg_facilities_column, agg_facilities_allcols, create_grid, read_pmeds, agg_pmeds_column,\
combine_pmeds_facilities, facilities_to_pmeds

def main():
    """ source activate pol """

    # combine_pmeds_facilities('Fresno', 'PM2.5T', damages=True, save=True)
    facilities_to_pmeds('LA', 'PM2.5T', save=True)
    exit()

if __name__ == "__main__":
    main()

def pop_dens(shpfn):
    pop = gpd.read_file('fresno_pop/fresnopop.shp').to_crs(epsg=4326)
    blocks = gpd.read_file(shpfn).to_crs(epsg=4326) # match crs of facilities
    pop['geometry'] = pop['geometry'].centroid # set all blocks to points
    # print(pop.iloc[0]['geometry'].centroid)
    # Gives table of facilities with column 'index_right' being which polygon each falls in
    facilities = gpd.sjoin(pop, blocks, op='within')
    name = 'POP10'
    # Sums pollutants for all points within a certain polygon
    total = facilities.groupby('index_right')[name].sum().to_frame() # need to convert Series to DF
    # Join counted values for each blockgroup back to blockgroups by blockgroup index and drop unmatched rows
    # Can also do blockgroups.merge(total, on="right_index")
    blocks = blocks.join(total).dropna()

    # pop['density'] = pop['POP10']/pop['geometry'].area
    # new_name = save_name('facilities/' + 'pop' + '_Fresno', '.shp')
    # blocks.to_file(new_name)
    return blocks

def pop_emissions(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    pop = pop_dens(shpfn)[['i', 'j', 'POP10']] # get rid of geometry for join
    em = calculate_total2(fn, shpfn, pol, name, points, scc_name, scc)

    # Add the two
    merged_df = em.merge(pop,  how='outer', left_on=['i','j'], right_on = ['i','j'])[[pol, 'POP10', 'geometry']] # can't hash on geometry
    # merged_df = gpd.GeoDataFrame(merged_df) # reconvert to gdf after merge
    merged_df['POP10'] = merged_df[pol] * merged_df['POP10']
    merged_df = merged_df[['i', 'j', 'POP10', 'geometry']]
    # new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    # merged_df.to_file(new_name)
    return merged_df

def pop_emissions2(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    df = pop_emissions('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name, scc)
    df2 = pop_emissions('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='Both', scc=[])


def norm_pop(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    pop = pop_dens(shpfn) # get rid of geometry for join
    max_v = max(pop['POP10'])
    pop['POP10'] = pop['POP10']/max_v
    em = calculate_total2(fn, shpfn, pol, name, points, scc_name, scc)[['i', 'j', pol]]

    # Add the two
    merged_df = pop.merge(em,  how='outer', left_on=['i','j'], right_on = ['i','j'])[[pol, 'POP10', 'geometry']] # can't hash on geometry
    # merged_df = gpd.GeoDataFrame(merged_df) # reconvert to gdf after merge
    merged_df['POP10'] = merged_df[pol] * merged_df['POP10']
    merged_df = merged_df[['POP10', 'geometry']]
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    merged_df.to_file(new_name)
    # return merged_df

"""
Folium Visualization
"""

# Creates choropleth folium map object
# http://andrewgaidus.com/leaflet_webmaps_python/
def add_choropleth(mapobj, gdf, id_field, value_field, units='', fill_color='YlOrRd', fill_opacity=0.6,
                    line_opacity=0.2, num_classes=5, classifier='Fisher_Jenks'):
    # 3 Pysal map classifiers to display data
    # Generate list of breakpoints using specified classification scheme. List of breakpoints will be input to choropleth function
    if classifier == 'Fisher_Jenks':
        threshold_scale=ps.esda.mapclassify.Fisher_Jenks(gdf[value_field], k = num_classes).bins.tolist()
    if classifier == 'Equal_Interval':
        threshold_scale=ps.esda.mapclassify.Equal_Interval(gdf[value_field], k = num_classes).bins.tolist()
    if classifier == 'Quantiles':
        threshold_scale=ps.esda.mapclassify.Quantiles(gdf[value_field], k = num_classes).bins.tolist()

    # Convert the GeoDataFrame to WGS84 coordinate reference system
    gdf_wgs84 = gdf.to_crs({'init': 'epsg:4326'})

    # Call Folium choropleth function, specifying the geometry as the WGS84 dataframe converted to GeoJSON, the data as
    # the GeoDataFrame, and the columns as the user-specified id and value fields.
    # key_on field refers to the id field within the GeoJSON string
    mapobj.choropleth(
                geo_data = gdf_wgs84.to_json(), # can also use folium.GeoJson()
                name = value_field,
                data = gdf,
                columns = [id_field, value_field],
                key_on = 'feature.properties.{}'.format(id_field),
                fill_color = fill_color,
                fill_opacity = fill_opacity,
                line_opacity = line_opacity,
                threshold_scale = threshold_scale,
                legend_name = value_field + units)
    return mapobj

# Add info popups for polygons
# https://github.com/python-visualization/folium/pull/376
def add_poly_info(mapobj, gdf, popup_field_list, name='Poly Info'):
    # Create a Folium feature group for this layer, since we will be displaying multiple layers
    poly_lyr = folium.FeatureGroup(name)
    j = json.loads(gdf.to_json())
    geojson = [{'type': j['type'], 'features': [f]} for f in j['features']]
    for gj in map(lambda gj: folium.GeoJson(gj, style_function=lambda feature: {
        # 'fillColor': feature['properties']['RGBA'],
        # 'color' : feature['properties']['RGBA'],
        'weight' : 1,
        'fillOpacity' : 0,
        }), geojson):
        label = '<br>'.join([field + ': ' + str(gj.data['features'][0]['properties'][field]) for field in popup_field_list])
        gj.add_child(folium.Popup(label)) # IFrame(label, width = 300, height = 100)
        gj.add_to(poly_lyr)
    # Add this point layer to the map object
    mapobj.add_child(poly_lyr)
    return mapobj

# Add points with clustering
def add_point_clusters(mapobj, gdf, popup_field_list):
    # Create empty lists to contain the point coordinates and the point pop-up information
    coords, popups = [], []
    # Loop through each record in the GeoDataFrame
    for i, row in gdf.iterrows():
        # Append lat and long coordinates to "coords" list
        coords.append([row.geometry.y, row.geometry.x])
        # Create a string of HTML code used in the IFrame popup
        # Join together the fields in "popup_field_list" with a linebreak between them
        label = '<br>'.join([field + ': ' + str(row[field]) for field in popup_field_list])
        # Append an IFrame that uses the HTML string to the "popups" list
        popups.append(IFrame(label, width = 300, height = 100))

    # Create a Folium feature group for this layer, since we will be displaying multiple layers
    pt_lyr = folium.FeatureGroup(name = 'Facilities')

    # Add the clustered points of crime locations and popups to this layer
    pt_lyr.add_child(MarkerCluster(locations = coords, popups = popups))

    # Add this point layer to the map object
    mapobj.add_child(pt_lyr)
    return mapobj

def add_top_point_clusters(mapobj, gdf, popup_field_list, n, col):
    gdf = gdf.sort_values(col, ascending = False).head(n)
    # Create empty lists to contain the point coordinates and the point pop-up information
    coords, popups = [], []
    # Loop through each record in the GeoDataFrame
    for i, row in gdf.iterrows():
        # Append lat and long coordinates to "coords" list
        coords.append([row.geometry.y, row.geometry.x])
        # Create a string of HTML code used in the IFrame popup
        # Join together the fields in "popup_field_list" with a linebreak between them
        label = '<br>'.join([field + ': ' + str(row[field]) for field in popup_field_list])
        # Append an IFrame that uses the HTML string to the "popups" list
        popups.append(IFrame(label, width = 300, height = 100))

    # Create a Folium feature group for this layer, since we will be displaying multiple layers
    pt_lyr = folium.FeatureGroup(name = 'Top ' + str(n) + ' Facilities')

    # Add the clustered points of crime locations and popups to this layer
    pt_lyr.add_child(MarkerCluster(locations = coords, popups = popups))

    # Add this point layer to the map object
    mapobj.add_child(pt_lyr)
    return mapobj

# Add points without clustering
def add_points(mapobj, gdf, popup_field_list):
    # Create a Folium feature group for this layer, since we will be displaying multiple layers
    pt_lyr = folium.FeatureGroup(name = 'pt_lyr')
    # Loop through each record in the GeoDataFrame
    for i, row in gdf.iterrows():
        label = '<br>'.join([field + ': ' + str(row[field]) for field in popup_field_list])
        pt_lyr.add_child(folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            popup=label,
            icon=folium.Icon(color='red', icon='info-sign')
        ))

    # Add this point layer to the map object
    mapobj.add_child(pt_lyr)
    return mapobj

def map_pollutant(name='PMT', basemap='OpenStreetMap'):
    # Create Fresno basemap specifying map center, zoom level, and using the default OpenStreetMap tiles
    # "Mapbox Bright" (Limited levels of zoom for free tiles)
    # "Mapbox Control Room" (Limited levels of zoom for free tiles)
    # "Stamen" (Terrain, Toner, and Watercolor)
    # "Cloudmade" (Must pass API key)
    # "Mapbox" (Must pass API key)
    # "CartoDB" (positron and dark_matter)
    m = folium.Map(location=[36.6, -119.6], zoom_start=9, tiles=basemap)

    # CES
    m = ces_layer(m)

    # Update basemap with choropleth
    pollutant_gdf = agg_points_col(name, 'fresno_industry_arb2.csv', 'ca_blockgroup/tl_2016_06_bg.shp')
    # make sure all values in column NAME have a value for viz
    pollutant_gdf[name].fillna(0)
    m = add_choropleth(m, pollutant_gdf, 'GEOID', name, units=' (Tons/yr)')
    m = add_poly_info(m, pollutant_gdf, [name], name + ' Info')

    # Update choropleth with facility clusters
    points = read_points('fresno_industry_arb2.csv')
    points['percent_tot'] = (points[name]/(points[name].sum()))*100
    m = add_point_clusters(m, points, ['FNAME', name, 'FSIC', 'percent_tot'])

    # Add top 40 polluting facilities
    m = add_top_point_clusters(m, points, ['FNAME', name, 'FSIC', 'percent_tot'], 40, name)

    # folium.plugins.MeasureControl().add_to(m)
    folium.LayerControl().add_to(m)

    # Save map
    m.save(save_name(name, ".html"))
    m

def most_polluting(pol, count=5):
    # Most-polluting FSIC code
    facilities = read_points('fresno_industry_arb2.csv')
    # facilities = facilities.groupby('FSIC')[pollutant].sum().to_frame()
    facilities = facilities.groupby('FSIC')[pol].agg({pol: 'sum', 'count': 'count'})
    print("Net Most-Polluting Facility SIC: \n{0}".format(facilities.loc[facilities[pol].idxmax()]))
    # Look up SIC: https://www.osha.gov/pls/imis/sic_manual.html
    return facilities.sort_values(by=pol, ascending=False).head(count)

def most_polluting_avg(pol, count=5):
    # Most-polluting FSIC code per "FSIC capita", so average pollution per facility type (FSIC code)
    facilities = read_points('fresno_industry_arb2.csv')
    # facilities = facilities.groupby('FSIC')[pollutant].sum().to_frame()
    facilities = facilities.groupby('FSIC')[pol].agg({pol: 'sum', 'count': 'count'})
    facilities['capita'] = facilities[pol]/facilities['count']
    print("Most-Polluting Facility Per Facility-Capita SIC: \n{0}".format(facilities.loc[facilities['capita'].idxmax()]))
    # Look up SIC: https://www.osha.gov/pls/imis/sic_manual.html
    return facilities.sort_values(by='capita', ascending=False).head(count)

def map_pop_dens():
    # https://www.census.gov/geo/maps-data/data/tiger-data.html
    pop = gpd.read_file('fresno_pop/fresnopop.shp')

    pop['density'] = pop['POP10']/pop['geometry'].area

    m = folium.Map(location=[36.6, -119.6], zoom_start=9, tiles='OpenStreetMap')

    # Update basemap with choropleth
    m = add_choropleth(m, pop, 'BLOCKID10', 'density', ' (population density)')

    folium.LayerControl().add_to(m)

    m.save(save_name('density', ".html"))
    m

def ces_layer(m):
    tracts = gpd.read_file('ces3shp/CES3Results.shp')
    # print("blocks shape before drop: {0}".format(blockgroups.shape))
    tracts = tracts[tracts['County'] == 'Fresno'] # Only keep blocks in that county
    # print(read_blocks('ces3shp/CES3Results.shp', '019').head())
    # print(tracts.columns)
    # print(tracts.head())
    # print(tracts.iloc[2]['Tract_1'])
    # m = folium.Map(location=[36.6, -119.6], zoom_start=9, tiles='OpenStreetMap')

    # Update basemap with choropleth
    m = add_choropleth(m, tracts, 'Tract_1', 'CIscore', ' (%)', fill_color='YlOrRd') # RdYlGn, BuGn, BuPu, GnBu, OrRd, PuBu, PuBuGn, PuRd, RdPu, YlGn, YlGnBu, YlOrBr, YlOrRd

    # folium.LayerControl().add_to(m)

    # m.save(save_name('ces', ".html"))
    return m



def map_grid(basemap='OpenStreetMap'):
    # gdf = gpd.read_file('grid.shp')
    gdf = gpd.read_file('MEDS/output/nox_Fresno_07d15col1.shp').to_crs(epsg=4326)
    # gdf_wgs84 = gdf.to_crs({'init': 'epsg:4326'})
    gj = json.loads(gdf.to_json())
    m = folium.Map(location=[36.6, -119.6], zoom_start=9, tiles=basemap)
    folium.GeoJson(gj, style_function=lambda feature: {
        # 'fillColor': feature['properties']['RGBA'],
        # 'color' : feature['properties']['RGBA'],
        'weight' : 1,
        'fillOpacity' : 0,
    }).add_to(m)

    # folium.plugins.MeasureControl().add_to(m)
    folium.LayerControl().add_to(m)

    # Save map
    m.save(save_name("name", ".html"))



def main():
    """ source activate pol """
    # agg_facilities_column('grid', 'Fresno', 'PM2.5T')
    # agg_pmeds_column('LA', 'PM2.5T', 'Both', 'percent', damages=True)

    # map_pollutant('PM2.5T', basemap='CartoDB Positron')
    # create_grid('grid')
    # map_grid()

    # print(most_polluting('PMT', 10))

    # print(most_polluting_avg('PMT', 10))

    # map_pop_dens()

    # convert_pmeds('MEDS/Fresno_07d15', 'grid.shp', 'pm', scc_name='LD', scc=[202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])
    # convert_pmeds('MEDS/Fresno_07d15', 'grid.shp', 'pm', scc_name='HD', scc=[302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,508,413,513,613,414,514,614,617,717,420,520,620])
    # normalize_pmeds('MEDS/Fresno_07d15', 'grid.shp', 'pm', scc_name='Both', scc=[])
    # normalize_total2('MEDS/Fresno_07d15', 'grid.shp', 'pm', scc_name='Fac')
    # calculate_total('MEDS/Fresno_07d15', 10, 'grid.shp', 'pm', scc_name='LD')
    # damages('MEDS/Fresno_07d15', 10, 'grid.shp', 'pm')
    # pop_dens('grid.shp')
    # agg_points_col2('PM2.5T', 10, 'grid.shp')
    # calculate_total2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='Both', scc=[])
    # damages2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='LD', scc=[202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])
    # damages3('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10)
    # damages2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='HD', scc=[302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,508,413,513,613,414,514,614,617,717,420,520,620])
    # damages2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='Both', scc=[])
    # normalize_total_damages('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='LD', scc=[202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])
    # normalize_total_damages('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='HD', scc=[302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,508,413,513,613,414,514,614,617,717,420,520,620])
    # normalize_total_damages2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='HD', scc=[])
    # pop_emissions('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='Both', scc=[])
    # norm_pop('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='Both', scc=[])
    # pop_emissions2('MEDS/Fresno_07d15', 'grid.shp', 'pm', 'PM2.5T', 10, scc_name='LD', scc=[202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])
    # grid = gpd.read_file('grid2.shp')
    # print(max(grid['i']))
    # print(max(grid['j']))
    # print(min(grid['i']))
    # print(min(grid['j']))
    # print(grid.head)
    # print(grid.shape)







# # facilities.intersects(blockgroups.unary_union)
#
# # Add column of points that are in polygons
# # https://stackoverflow.com/questions/27606924/count-number-of-points-in-multipolygon-shapefile-using-python/27608968#27608968
# def add_intersects():
#     blockgroups['PMT'] = 0
#     for i, poly in blockgroups.iterrows():
#         # Keep a list of points in this poly
# #         pts_in_this_poly = []
#         sum_pts = 0
#
#         # Now loop over all points with index j.
#         for j, pt in facilities.iterrows():
#             if poly.geometry.contains(pt.geometry):
#                 # Then it's a hit! Add it to the list,
#                 # and drop it so we have less hunting.
# #                 pts_in_this_poly.append(pt.geometry)
#                 sum_pts += pt['PMT']
#                 pts = facilities.drop([j])
#
#         # We could do all sorts, like grab a property of the
#         # points, but let's just append the number of them.
# #         pts_in_polys.append(len(pts_in_this_poly))
#         blockgroups['PMT'][i] = sum_pts
#
# # add_intersects()
#
# # print(blockgroups.head())



# # # blockgroups['PMT'][0] = 1
# # blockgroups.set_value('PMT', 0, 10)
# # set(blockgroups['PMT'])
# blockgroups = add_pollutants()
# blockgroups.columns
# # print(blockgroups.iloc[0])
# print(blockgroups.loc['060190001001'])
#
# blockgroups = add_pollutant('SOXT')
# print(blockgroups.loc[blockgroups['GEOID'] == '060190001001'])



# TODO:
# print(set(df['FSIC']))
# for i in range(len(df)):
#     if df['FSIC'][i] > 2899 and df['FSIC'][i] < 3000:
#         print(df.iloc[[i]])
# find most polluting fsic per block group
# just do nox, sox, pm10
# just dairy

# print(df.head())
# map



# for i in range(len(df)):
#     map.draw({'x': df['lon'][i], 'y': df['lat'][i]})
