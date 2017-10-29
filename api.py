"""
Imports
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

"""
Helpers
"""

# Pick a nonconflicting save name for a file to not overwrite previous ones.
def save_name(name, extension):
    i = 0
    fname = name
    while os.path.isfile(fname + extension):
        fname = name + str(i)
        i += 1
    return fname + extension

# Reads blocks into a geopandas dataframe
def read_blocks(name, county_code):
    blockgroups = gpd.read_file(name)
    # print("blocks shape before drop: {0}".format(blockgroups.shape))
    blockgroups = blockgroups[blockgroups['COUNTYFP'] == county_code] # Only keep blocks in that county
    # print("blocks shape after drop: {0}".format(blockgroups.shape))
    return blockgroups

# Converts csv of lats/lons to a geopandas dataframe
# https://gis.stackexchange.com/questions/114066/handling-kml-csv-with-geopandas-drivererror-unsupported-driver-ucsv
def read_points(name):
    facilities = pd.read_csv(name)
    # print("facilities shape before drop: {0}".format(facilities.shape))
    facilities = facilities[pd.notnull(facilities['lat'])] # drop ungeocoded rows; can also use np.isfinite
    facilities['geometry'] = facilities.apply(lambda z: Point(z.lon, z.lat), axis=1)
    facilities = gpd.GeoDataFrame(facilities)
    facilities = facilities[facilities.within(read_blocks('ca_blockgroup/tl_2016_06_bg.shp', '019').unary_union)] # drop rows outside region
    # facilities = gpd.GeoDataFrame(facilities)
    # print("facilities shape after drop: {0}".format(facilities.shape)) # one more column since we added geometry
    return facilities

# Adds up column NAME for all POINTS in each POLYGON
def agg_points_col(name, points, polygons):
    blockgroups = read_blocks(polygons, '019')
    # print("blocks shape before add: {0}".format(blockgroups.shape))
    facilities = read_points(points)
    # make sure all values in column NAME have a value in case
    facilities[name].fillna(0)
    facilities.crs = blockgroups.crs # make facilities projection that of blockgroups
    # Gives table of facilities with column 'index_right' being which polygon each falls in
    facilities = gpd.sjoin(facilities, blockgroups, op='within')
    # Sums pollutants for all points within a certain polygon
    total = facilities.groupby('index_right')[name].sum().to_frame() # need to convert Series to DF
    # Join counted values for each blockgroup back to blockgroups by blockgroup index and drop unmatched rows
    # Can also do blockgroups.merge(total, on="right_index")
    blockgroups = blockgroups.join(total).dropna()
    # print("blocks shape after add: {0}".format(blockgroups.shape))
    return blockgroups

# Adds up all columns from all points in each blockgroup
def agg_points_allcols(points, polygons):
    blockgroups = read_blocks(polygons, '019')
    # print("blocks shape before add: {0}".format(blockgroups.shape))
    facilities = read_points(points)
    facilities.crs = blockgroups.crs # make facilities projection that of blockgroups
    blockgroups = gpd.sjoin(blockgroups, facilities)
    blockgroups = blockgroups.dissolve(by='GEOID', aggfunc='sum')
    # Alternative is blockgroups.groupby('geometry')['PMT'].sum()
    # print("blocks shape after add: {0}".format(blockgroups.shape))
    return blockgroups

# Narrows blocks from a state to a specified county and saves if out_file is not None
def narrow_state_blocks(county_code, in_file, out_file=None, extension=".shp"):
    blocks = gpd.read_file(in_file + extension)
    blocks = blocks[blocks['COUNTYFP10'] == county_code] # Only keep blocks in that county
    if out_file:
        blocks.to_file(save_name(out_file, extension))
    return blocks

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

"""
Main Visualization Functions
"""

def map_pollutant(name='PMT', basemap='OpenStreetMap'):
    # Create Fresno basemap specifying map center, zoom level, and using the default OpenStreetMap tiles
    # “Mapbox Bright” (Limited levels of zoom for free tiles)
    # “Mapbox Control Room” (Limited levels of zoom for free tiles)
    # “Stamen” (Terrain, Toner, and Watercolor)
    # “Cloudmade” (Must pass API key)
    # “Mapbox” (Must pass API key)
    # “CartoDB” (positron and dark_matter)
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

# https://gis.stackexchange.com/questions/54119/creating-square-grid-polygon-shapefile-with-python?rq=1
# measurements in meters
def create_grid(name="test.shp", xmin=-684000, xmax=600000, ymin=-564000, ymax=600000, gridHeight=4000, gridWidth=4000):
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = math.ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = math.ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(name):
        os.remove(name)
    outDataSource = outDriver.CreateDataSource(name)
    outLayer = outDataSource.CreateLayer(name,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature.Destroy

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Close DataSources
    outDataSource.Destroy()

def create_grid2(name="test", minx=-684000, maxx=600000, miny=-564000, maxy=600000, dx=4000, dy=4000):

    nx = int(math.ceil(abs(maxx - minx)/dx))
    ny = int(math.ceil(abs(maxy - miny)/dy))

    w = shp.Writer(shp.POLYGON)
    w.autoBalance = 1
    w.field("ID")
    w.field('i')
    w.field('j')
    id=0

    for i in range(ny):
        for j in range(nx):
            id+=1
            vertices = []
            parts = []
            vertices.append([min(minx+dx*j,maxx),max(maxy-dy*i,miny)])
            vertices.append([min(minx+dx*(j+1),maxx),max(maxy-dy*i,miny)])
            vertices.append([min(minx+dx*(j+1),maxx),max(maxy-dy*(i+1),miny)])
            vertices.append([min(minx+dx*j,maxx),max(maxy-dy*(i+1),miny)])
            parts.append(vertices)
            w.poly(parts)
            w.record(id,i,j)

    w.save(name)
    # projection info
    # prj = open("%s.prj" % name, "w")
    # epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6370000,6370000]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    # epsg = 'GEOGCS["LCC",DATUM["LCC",SPHEROID["LCC",6370000,6370000]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
    # prj.write(epsg)
    # prj.close()

def map_grid(basemap='OpenStreetMap'):
    gdf = gpd.read_file('test.shp')
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

def convert_pmeds(fn, shpfn):
    # Read in fixed-width columns
    df = pd.read_fwf(fn, widths=[8, 14, 14, 3, 3, 9, 5, 2, 5, 2, 2, 3, 3, 5], skiprows=0, header=None)
    # Read comma-separated columns; assumes no commas in first 79 characters
    df2 = pd.read_csv(fn, header=None)
    df2[df2.columns[0]] = df2[df2.columns[0]].map(lambda x: x[78:])
    # Join dataframes
    df = pd.concat([df, df2], axis=1)
    # Name column
    df.columns = ['scenario', 'sic', 'scc', 'i', 'j', 'facid', 'stackid', 'county',
    'julianday', 'beginhour', 'endhour', 'airbasin', 'subcounty', 'elevation', 'co', 'nox', 'sox', 'tog', 'pm', 'nh3']
    # print(df.head())
    # Add up individual facilities' pm
    df['pm'] = df.groupby(['i', 'j'])['pm'].transform('sum')
    df = df.drop_duplicates(['i', 'j'])

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    gdf = gpd.read_file(shpfn)
    # print(gdf.shape)

    merged_df = gdf.merge(df,  how='left', left_on=['j','i'], right_on = ['i','j'])[['pm', 'geometry']] # columns to keep
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)


def main():
    """ source activate pol """
    # map_pollutant('PM2.5T', basemap='CartoDB Positron')
    # create_grid2('test2')
    # map_grid()

    # print(most_polluting('PMT', 10))

    # print(most_polluting_avg('PMT', 10))

    # map_pop_dens()

    convert_pmeds('MEDS/Fresno_07d31', 'test2.shp')


if __name__ == "__main__":
    main()




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
