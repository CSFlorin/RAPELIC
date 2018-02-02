"""
California Industrial and Auto Pollution API Helpers
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


"""
Helpers
"""

def save_name(name, extension=''):
    """ Return nonconflicting save name including extension for a file """
    i = 0
    fname = name
    while os.path.isfile(fname + extension):
        fname = name + '_' + str(i)
        i += 1
    return fname + extension

def county_code(county, file_type):
    """
    Return the properly formatted county code string for a county
        county: 'LA', 'Fresno'
        file_type: 'census_area', 'arb_facilities'
    """
    if file_type == 'census_area':
        # All codes: https://www2.census.gov/geo/docs/reference/codes/files/st06_ca_cou.txt
        if county == 'LA':
            return '037'
        if county == 'Fresno':
            return '019'
    elif file_type == 'arb_facilities':
        if county == 'LA':
            return '19'
        if county == 'Fresno':
            return '10'

def read_census_shape(area, county=''):
    """
    Read California census area shapefile into a GeoDataFrame and return it filtered to a county if specified
        area: 'blockgroup', 'block'
        county: 'LA', 'Fresno'
    """
    if area == 'blockgroup':
        area = gpd.read_file('data/ca_blockgroup/tl_2016_06_bg.shp')
    elif area == 'block':
        area = gpd.read_file('data/ca_block/tl_2016_06_tabblock10.shp')
    else:
        raise ValueError('area must be LA or Fresno')

    if county: # Only keep area in that county if specified
        area = area[area['COUNTYFP'] == county_code(county, 'census_area')]
    return area

def read_facilities(county):
    """
    Read csv of lats/lons and return a GeoDataFrame
        county: 'LA', 'Fresno'
    """
    facilities = pd.read_csv('data/ARB_OUT/' + county_code(county, 'arb_facilities') + '.csv')
    # https://gis.stackexchange.com/questions/114066/handling-kml-csv-with-geopandas-drivererror-unsupported-driver-ucsv
    facilities['geometry'] = facilities.apply(lambda z: Point(z.lon, z.lat), axis=1)
    facilities = gpd.GeoDataFrame(facilities)
    facilities = facilities[facilities.within(read_census_shape('blockgroup', county).unary_union)] # Drop rows outside region in case any were misgeocoded
    return facilities

def agg_facilities_column(area, county, column, output='percent', damages=True, save=False):
    """
    Return geopandas dataframe with values of column added up for all facilities in each polygon
        area: 'blockgroup', 'block', 'grid'
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PM10T', 'PMT', 'SOXT', 'NOXT', 'COT', 'ROGT', 'TOGT'
        output: 'percent' for percent that each cell is of all cells, 'actual'
        damages: If output == 'percent', whether to calculate for cell damages or just column values
        save: Whether to save frame to shapefile
    """
    if output == 'actual' and damages:
        raise ValueError('Cannot calculate damages if outputting actual column values')

    if area == 'grid':
        area = gpd.read_file('data/MEDS/grid.shp').to_crs(epsg=4326) # Match crs of facilities
    else:
        area = read_census_shape(area, county)

    facilities = read_facilities(county)
    # Convert facilities crs to match area's
    facilities.crs = area.crs
    # facilities.crs = {'init': 'epsg:4326', 'no_defs': True} # crs for naive geometries per https://github.com/geopandas/geopandas/issues/245
    # Ensure all values in column have a value for addition
    facilities[column].fillna(0)
    # Give table of facilities with new column 'index_right' being which polygon each falls in
    facilities = gpd.sjoin(facilities, area, op='within')
    # Sum pollutants for all points within each polygon
    total = facilities.groupby('index_right')[column].sum().to_frame() # Convert Series to DataFrame
    # Join counted values for each area back to area by index and drop unmatched rows
    area = area.join(total).dropna() # Can also do blockgroups.merge(total, on="right_index")

    if output != 'actual':
        sum_v = sum(area[column])
        if damages:
            if county == 'Fresno':
                area[column] = (area[column]/sum_v)*19900000 # Total PM2.5 damages for Fresno
            elif county == 'LA':
                area[column] = (area[column]/sum_v)*243000000 # Total PM2.5 damages for LA
        else:
            # Divide each cell by total
            area[column] = (area[column]/sum_v)*100

    if save:
        new_name = save_name('out/facilities/' + county + column, '.shp')
        area.to_file(new_name)

    return area

def agg_facilities_allcols(area, county):
    """
    Return geopandas dataframe with values of each column added up for all facilities in each polygon
    No filtering down to just one column
        area: 'blockgroup', 'block'
        county: 'LA', 'Fresno'
    """
    area = read_census_shape(area, county)
    facilities = read_facilities(county)
    # Convert facilities crs to match area's
    facilities.crs = area.crs
    area = gpd.sjoin(area, facilities)
    area = area.dissolve(by='GEOID', aggfunc='sum') # Alternative is blockgroups.groupby('geometry')['PMT'].sum()
    return area

def create_grid(name='grid', minx=-684000, maxx=600000, miny=-564000, maxy=600000, dx=4000, dy=4000):
    """ Write PMEDS grid shapefile composed of squares, where i and j correspond to PMEDS """
    nx = int(math.ceil(abs(maxx - minx)/dx))
    ny = int(math.ceil(abs(maxy - miny)/dy))
    shift = (maxy - miny)/dy
    w = shp.Writer(shp.POLYGON)
    w.autoBalance = 1
    w.field('ID', 'C')
    w.field('i', 'N') # Make int instead of object
    w.field('j', 'N')
    square_id = 0
    for j in range(0, ny): # 0-index in loop for calculations
        for i in range(0, nx):
            vertices = []
            parts = []
            vertices.append([min(minx+dx*i,maxx), max(maxy-dy*j,miny)])
            vertices.append([min(minx+dx*(i+1),maxx), max(maxy-dy*j,miny)])
            vertices.append([min(minx+dx*(i+1),maxx), max(maxy-dy*(j+1),miny)])
            vertices.append([min(minx+dx*i,maxx), max(maxy-dy*(j+1),miny)])
            parts.append(vertices)
            w.poly(parts)
            w.record(square_id,i+1,shift-j) # 1-index when storing i and j
            square_id += 1
    w.save(save_name(name))

def read_pmeds(county):
    """
    Return dataframe of PMEDS file
        county: 'LA', 'Fresno'
    """
    fn = 'data/MEDS/' + county + '_07d15'
    # Read in fixed-width columns
    df = pd.read_fwf(fn, widths=[8, 14, 14, 3, 3, 9, 5, 2, 5, 2, 2, 3, 3, 5], skiprows=0, header=None)
    # Read comma-separated columns; assumes no commas in first 79 characters since using read_csv
    df2 = pd.read_csv(fn, header=None)
    df2[df2.columns[0]] = df2[df2.columns[0]].map(lambda x: x[78:])
    # Concatenate dataframes
    df = pd.concat([df, df2], axis=1)
    # Name columns
    df.columns = ['scenario', 'sic', 'SCC', 'i', 'j', 'facid', 'stackid', 'county',
    'julianday', 'beginhour', 'endhour', 'airbasin', 'subcounty', 'elevation', 'COT',
    'NOXT', 'SOXT', 'TOGT', 'PMT', 'NH3']

    # Calculate pm2.5 for pmeds
    ratio_df = pd.read_csv('data/MEDS/rogpm.11jan2017.txt')
    ratio_df = ratio_df.loc[ratio_df['Year'] == 2010] # year is int, and match with year of pmeds
    # Match scc, year in ratio_df to scc in df, scc is int, sic is float
    df = df.merge(ratio_df,  how='inner', on=['SCC']) # columns to keep
    df['PM2.5T'] = df['PMT'] * df['PM2_5/TPM'] # Multiply pm by conversion factor for pm2.5
    # Fill NaN values with 0
    df.fillna(0, inplace=True)
    return df

def agg_pmeds_column(county, column, scc='', output='percent', updated=True, damages=False, save=False):
    """
    Return geodataframe of PMEDS column added up for each grid cell
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'NH3', 'TOGT'
        scc: 'LD', 'HD', '' for all
        output: 'percent', 'actual'
        updated: If output == 'actual', True to update to 2015 values
        damages: If output == 'percent', True to calculate cell damages and False for just cell percents
        save: Whether to save frame to shapefile
    """
    if output == 'actual' and damages:
        raise ValueError('Cannot calculate damages if outputting actual column values')

    df = read_pmeds(county)

    # Only keep SCCs we want, if specified
    if scc == 'LD':
        df = df[df['SCC'].isin([202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])]
    elif scc == 'HD':
        df = df[df['SCC'].isin([302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,508,413,513,613,414,514,614,617,717,420,520,620])]

    # Add up individual squares' pm
    df[column] = df.groupby(['i', 'j'])[column].transform('sum') # Add each row containing same i, j, but keep dupes
    df = df.drop_duplicates(['i', 'j']) # Get rid of dupes

    if output == 'actual':
        if updated:
            # Divide each cell by total
            sum_v = sum(df[column])
            df[column] = df[column]/sum_v
            # Update to 2015 vals using tons/day
            if county == 'Fresno':
                # https://www.arb.ca.gov/app/emsinv/2017/emssumcat_query.php?F_YR=2015&F_DIV=-4&F_SEASON=A&SP=SIP105ADJ&F_AREA=CO&F_CO=10&F_COAB=#7
                if scc == 'LD':
                    df[column] = df[column]*.44*365
                elif scc == 'HD':
                    df[column] = df[column]*.65*365
                else:
                    df[column] = df[column]*1.09*365
            elif county == 'LA':
                # https://www.arb.ca.gov/app/emsinv/2017/emssumcat_query.php?F_YR=2015&F_DIV=-4&F_SEASON=A&SP=SIP105ADJ&F_AREA=CO&F_CO=19&F_COAB=#7
                if scc == 'LD':
                    df[column] = df[column]*4.9*365
                elif scc == 'HD':
                    df[column] = df[column]*2.05*365
                else:
                    df[column] = df[column]*6.95*365
        else:
            df[column] = df[column]*0.00110231*24*365 # Convert kg/hr to tons/yr after aggregation
    else:
        if damages:
            sum_v = sum(df[column])
            if county == 'Fresno':
                if scc == 'LD':
                    df[column] = (df[column]/sum_v)*9741393 # LD damages in dollars
                elif scc == 'HD':
                    df[column] = (df[column]/sum_v)*18369483
                else:
                    df[column] = (df[column]/sum_v)*28110876
            elif county == 'LA':
                if scc == 'LD':
                    df[column] = (df[column]/sum_v)*457609951
                elif scc == 'HD':
                    df[column] = (df[column]/sum_v)*281522511
                else:
                    df[column] = (df[column]/sum_v)*739132462
        else:
            # Divide each cell by total
            sum_v = sum(df[column])
            df[column] = (df[column]/sum_v)*100

    # Join with grid shapefile
    gdf = gpd.read_file('data/MEDS/grid.shp') # gdf shape 93411 by 4

    # Keep only right or inner
    merged_df = gdf.merge(df,  how='right', on=['i','j'])[[column, 'i', 'j', 'geometry', 'ID']] # Columns to keep, shape 599 by 2
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())

    if save:
        new_name = save_name('out/MEDS/' + county + column + scc, '.shp')
        merged_df.to_file(new_name)
        # Add projection file
        from shutil import copy2
        copy2('data/MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

    return merged_df

def combine_pmeds_facilities(county, column, scc='', output='percent', updated=True, damages=False, save=False):
    """
    Return geodataframe of column added from PMEDS and facilities for each grid cell
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'TOGT'
        scc: 'LD', 'HD', '' for all
        output: 'percent', 'actual'
        updated: If output == 'actual', True to update to 2015 values
        damages: If output == 'percent', True to calculate cell damages and False for just cell percents
        save: Whether to save frame to shapefile
    """
    if output == 'actual' and damages:
        raise ValueError('Cannot calculate damages if outputting actual column values')

    # Calculate raw PMEDS frame if not damages, else damages frame
    pmeds = agg_pmeds_column(county, column, scc, 'actual' if not damages else 'percent', updated, damages)
    # Calculate raw facilities frame if not damages, else damages frame
    facilities = agg_facilities_column('grid', county, column, 'actual' if not damages else 'percent', damages)
    df = pd.concat([pmeds, facilities])
    # print(df.sort_values(['i', 'j'], ascending=[True, True]).head())
    # Add up individual squares' pm
    df[column] = df.groupby(['i', 'j'])[column].transform('sum') # Add each row containing same i, j, but keep dupes
    df = df.drop_duplicates(['i', 'j']) # Get rid of dupes

    if output == 'percent' and not damages:
        # Divide each cell by total
        sum_v = sum(df[column])
        df[column] = (df[column]/sum_v)*100

    if save:
        new_name = save_name('out/combined/' + county + column + scc, '.shp')
        df.to_file(new_name)
        # Add projection file
        from shutil import copy2
        copy2('data/MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

    return df
