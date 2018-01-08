"""
California Industrial and Auto Pollution API
Author: Florin Langer, UC Berkeley Renewable and Appropriate Energy Lab with Lawrence Berkeley National Laboratory
Sources: California Air Resources Board
"""


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

def agg_facilities_column(area, county, column, output='percent', damages=True):
    """
    Return geopandas dataframe with values of column added up for all facilities in each polygon
        area: 'blockgroup', 'block', 'grid'
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PM10T', 'PMT', 'SOXT', 'NOXT', 'COT', 'ROGT', 'TOGT'
        output: 'percent', 'actual'
        damages: If output == 'percent', whether to calculate cell damages or just cell percents
    """
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
    # area.crs = facilities.crs # make facilities projection that of blockgroups
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

    new_name = save_name('out/facilities/' + county + column, '.shp')
    area.to_file(new_name)
    return area

def agg_facilities_allcols(area, county):
    """
    Return geopandas dataframe with values of each column added up for all facilities in each polygon
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
    """ Write PMEDS grid shapefile composed of squares """
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

def agg_pmeds_column(county, column, scc='', output='percent', updated=True, damages=False):
    """
    Return geodataframe of PMEDS column added up for each grid cell
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'NH3', 'TOGT'
        scc: 'LD', 'HD', ''
        output: 'percent', 'actual'
        updated: If output == 'actual', True to update to 2015 values
        damages: If output == 'percent', True to calculate cell damages and False for just cell percents
    """
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
                    df[column] = (df[column]/sum_v)*9741393 # LD damages
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
    new_name = save_name('out/MEDS/' + county + column + scc, '.shp')
    merged_df.to_file(new_name)

    # Add projection file
    from shutil import copy2
    copy2('data/MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))
    return merged_df

# TODO: Fix following code

""" Return dataframe with total pm in a column """
def add_pmeds_facilities(fn, coid, shpfn, pol):
    pmeds = read_pmeds(fn)
    # Add up individual squares' pm
    pmeds[pol] = pmeds.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    pmeds = pmeds.drop_duplicates(['i', 'j']) # get rid of dupe rows
    # Normalize PMEDS from 0 to 1
    pol_max = max(pmeds[pol])
    pmeds[pol] = pmeds[pol]/pol_max
    # Adjust PMEDS for more accurate data and to tons/year
    pmeds[pol] = pmeds[pol].map(lambda x: 365*1.09*x)

    grid = gpd.read_file(shpfn).to_crs(epsg=4326) # match crs of facilities
    facilities = read_points('ARB_OUT/' + str(coid) + '.csv')
    facilities.crs = {'init': 'epsg:4326', 'no_defs': True} # the crs for naive geometries per https://github.com/geopandas/geopandas/issues/245

    pol2 = pol.upper() + '2.5T' if pol == 'pm' else pol.upper() + 'T' # pol name in facilities
    facilities[pol2].fillna(0) # make sure all values in column NAME have a value in case
    grid.crs = facilities.crs # make facilities projection that of blockgroups
    # Gives table of facilities with column 'index_right' being which polygon each falls in
    facilities = gpd.sjoin(facilities, grid, how='right', op='within')
    # Sums pol2 for all facilities within a certain square
    facilities = facilities.groupby('index_right')[pol2].sum().to_frame() # need to convert Series to DF, this is just index_right, PM2.5T
    # Join counted values for each blockgroup back to blockgroups by blockgroup index and drop unmatched rows
    # Can also do blockgroups.merge(total, on="right_index")
    facilities = grid.join(facilities) # join back with grid for i, j values and geometry
    gdf = facilities.merge(pmeds, how='outer', on=['i','j'])[['i', 'j', 'geometry', pol2, pol]]
    gdf[pol] = gdf.apply(lambda row: row[pol] + row[pol2], axis=1)
    gdf.to_file('test2.shp')
    gdf = gdf[pd.notnull(gdf[pol])]

    return gdf[['i', 'j', 'geometry', pol]]

""" Get percent that SCC_NAME (HD, LD, or I) is of total of POL
    COID is 19 for LA, 10 for Fresno
"""
def calculate_total(fn, coid, shpfn, pol, scc_name=''):
    df = read_pmeds(fn)

    # Codes from PMEDS SCC Emission Categories Chart
    if scc_name == 'HD':
        numerator = df[df['scc'].isin([302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,508,413,513,613,414,514,614,617,717,420,520,620])].copy()
    else: # LD
        numerator = df[df['scc'].isin([202,203,204,205,206,207,209,210,211,212,213,214,215,216,218,219,221,808,813,814,817,820])].copy()

    # Add up individual squares' pm
    df[pol] = df.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    df = df.drop_duplicates(['i', 'j']) # get rid of dupe rows
    numerator[pol] = numerator.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    numerator = numerator.drop_duplicates(['i', 'j']) # get rid of dupes

    # Normalize scc_name over total pol for pmeds
    numerator = numerator.merge(df,  how='inner', on=['i','j'])[['i', 'j', pol+'_x', pol+'_y']]
    df_max = max(df[pol])
    numerator[pol] = numerator.apply(lambda row: row[pol+'_x']/df_max, axis=1) # row[pol+'_y']
    # Multiply by appropriate mobile sources from https://www.arb.ca.gov/app/emsinv/2017/emssumcat_query.php?F_YR=2015&F_DIV=-4&F_SEASON=A&SP=SIP105ADJ&F_AREA=CO&F_CO=10&F_COAB=#7
    if scc_name == 'LD':
        numerator[pol] = numerator[pol].map(lambda x: 365*.35*x) # and convert from tons/day to tons/year
    else: # HD
        numerator[pol] = numerator[pol].map(lambda x: 365*.74*x)
    numerator = numerator[[pol, 'i', 'j']]

    total = add_pmeds_facilities(fn, coid, shpfn, pol)
    # Keep only right or inner
    merged_df = total.merge(numerator,  how='inner', on=['i','j'])[[pol + '_x', pol + '_y', 'geometry']] # columns to keep, shape 599 by 2
    merged_df[pol] = merged_df.apply(lambda row: row[pol+'_y']/row[pol+'_x'], axis=1) # numerator/total
    # merged_df = merged_df[pd.notnull(merged_df[pol])]
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name + '_total', '.shp')
    merged_df[[pol, 'geometry']].to_file(new_name)

    # No need to add projection file since we reprojected in geopandas
    # Add projection file
    # from shutil import copy2
    # copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

""" Total emissions """
def calculate_total2(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    """ PMEDS """
    df = read_pmeds(fn)
    # only keep SCCs we want, if specified
    if scc:
        print(df.shape)
        df = df[df['scc'].isin(scc)]
        print(df.shape)

    # Add up individual sources' pm
    df[pol] = df.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    df = df.drop_duplicates(['i', 'j']) # get rid of dupes

    # Get max total square
    df2 = read_pmeds(fn)
    df2[pol] = df2.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    df2 = df2.drop_duplicates(['i', 'j']) # get rid of dupes
    max_v = max(df2[pol])

    df[pol] = df[pol]/max_v

    # Update to 2015 values
    if scc_name == 'LD':
        df[pol] = df[pol]*.44*365
    elif scc_name == 'HD':
        df[pol] = df[pol]*.65*365
    elif scc_name == 'Both':
        df[pol] = df[pol]*1.09*365

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    # # gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # print(gdf.dtypes)

    # NOTE: Make sure dtypes of the columns in both tables are the same
    # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # x/j is 1 to 321, y/i is 1 to 291
    # Keep only right or inner
    # # merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, 'geometry']] # columns to keep, shape 599 by 2
    # now have raw values for pmeds
    # new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    # merged_df.to_file(new_name)
    # exit()
    """ FACILITIES """
    blockgroups = gpd.read_file(shpfn).to_crs(epsg=4326) # match crs of facilities
    # print("blocks shape before add: {0}".format(blockgroups.shape))
    facilities = read_points('ARB_OUT/' + str(points) + '.csv')
    facilities.crs = {'init': 'epsg:4326', 'no_defs': True} # the crs for naive geometries per https://github.com/geopandas/geopandas/issues/245

    # make sure all values in column NAME have a value in case
    facilities[name].fillna(0)
    blockgroups.crs = facilities.crs # make facilities projection that of blockgroups
    # Gives table of facilities with column 'index_right' being which polygon each falls in
    facilities = gpd.sjoin(facilities, blockgroups, op='within')
    # Sums pollutants for all points within a certain polygon
    total = facilities.groupby('index_right')[name].sum().to_frame() # need to convert Series to DF
    # Join counted values for each blockgroup back to blockgroups by blockgroup index and drop unmatched rows
    blockgroups = blockgroups.join(total).dropna()
    # now have raw values for facilities
    # print(blockgroups.dtypes)
    # blockgroups = blockgroups[['i', 'j', name]] # get rid of geometry col

    print(df.loc[(df['i'] == 165) & (df['j'] == 138)])
    print(blockgroups.loc[(blockgroups['i'] == 165) & (blockgroups['j'] == 138)])
    # exit()
    # Add the two
    # merged_df = merged_df.merge(blockgroups,  how='outer', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, name, 'geometry']] # can't hash on geometry
    merged_df = blockgroups.merge(df,  how='outer', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, name]] # can't hash on geometry
    print(merged_df.loc[(merged_df['i'] == 165) & (merged_df['j'] == 138)])
    # exit()
    # merged_df = gpd.GeoDataFrame(merged_df) # reconvert to gdf after merge
    # merged_df[pol] = merged_df[pol] + merged_df[name]
    merged_df.fillna(0, inplace=True) # when adding nan to int, result was nan
    merged_df[pol] = merged_df.apply(lambda row: row[pol] + row[name], axis=1)
    print(merged_df.loc[(merged_df['i'] == 165) & (merged_df['j'] == 138)])
    # exit()
    # join with grid
    gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    merged_df = gdf.merge(merged_df,  how='inner', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, 'geometry']]
    # merged_df = merged_df[['i', 'j', pol, 'geometry']]

    return merged_df

    # new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    # merged_df.to_file(new_name)
    #
    # # Add projection file
    # from shutil import copy2
    # copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

""" normalized ld and hd """
def normalize_total(fn, shpfn, pol, scc_name='', scc=[]):
    df = new_pmeds(fn, shpfn, pol, scc_name, scc)

    # Get max total square
    df2 = calculate_total2(fn, shpfn, pol, 'PM2.5T', 10, scc_name='Both', scc=[])
    max_v = max(df2[pol])

    df[pol] = df[pol]/max_v

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    # gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # # print(gdf.dtypes)
    #
    # # NOTE: Make sure dtypes of the columns in both tables are the same
    # # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # # x/j is 1 to 321, y/i is 1 to 291
    # # Keep only right or inner
    # merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[[pol, 'geometry']] # columns to keep, shape 599 by 2
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    df.to_file(new_name)

    # Add projection file
    from shutil import copy2
    copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

""" for finding normalized facilities """
def normalize_total2(fn, shpfn, pol, name='PM2.5T', scc_name=''):
    df = agg_points_col2(name, 10, shpfn)

    # Get max total square
    df2 = calculate_total2(fn, shpfn, pol, name, 10, scc_name='Both', scc=[])
    max_v = max(df2[pol])

    df[name] = df[name]/max_v

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    # gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # # print(gdf.dtypes)
    #
    # # NOTE: Make sure dtypes of the columns in both tables are the same
    # # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # # x/j is 1 to 321, y/i is 1 to 291
    # # Keep only right or inner
    # merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[[pol, 'geometry']] # columns to keep, shape 599 by 2
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    df.to_file(new_name)

    # Add projection file
    # from shutil import copy2
    # copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))



""" PMEDS damages """
def damages2(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    """ PMEDS """
    df = read_pmeds(fn)
    # only keep SCCs we want, if specified
    if scc:
        print(df.shape)
        df = df[df['scc'].isin(scc)]
        print(df.shape)

    # Add up individual sources' pm
    df[pol] = df.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    df = df.drop_duplicates(['i', 'j']) # get rid of dupes

    # Get total square
    df2 = read_pmeds(fn)
    if scc:
        df2 = df2[df2['scc'].isin(scc)]
    df2[pol] = df2.groupby(['i', 'j'])[pol].transform('sum') # adds each row containing same i, j, but keeps dupes
    df2 = df2.drop_duplicates(['i', 'j']) # get rid of dupes

    # df = df.merge(df2,  how='inner', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', 'pm_x', 'pm_y']]
    df2.fillna(0, inplace=True)
    max_v = sum(df2[pol])
    print(max_v)
    df[pol] = df[pol]/max_v

    # Caclulate damages
    if scc_name == 'LD':
        df[pol] = df[pol]*(1572107734.352)
    elif scc_name == 'HD':
        df[pol] = df[pol]*(2298989169-(1572107734.352))
    elif scc_name == 'Both':
        df[pol] = df[pol]*(2298989169)

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # print(gdf.dtypes)

    # NOTE: Make sure dtypes of the columns in both tables are the same
    # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # x/j is 1 to 321, y/i is 1 to 291
    # Keep only right or inner
    merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, 'geometry']] # columns to keep, shape 599 by 2
    # now have raw values for pmeds
    return merged_df
    # n = int(merged_df.shape[0] * .25)
    # merged_df = merged_df.nlargest(n, pol)
    # new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    # merged_df.to_file(new_name)
    # exit()

""" all damages """
def damages3(fn, shpfn, pol, name, points, scc_name='Both', scc=[]):
    fac = damages(fn, points, shpfn, pol)
    pmeds = damages2(fn, shpfn, pol, name, points, scc_name, scc)
    # add
    merged_df = fac.merge(pmeds,  how='outer', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', 'PM2.5 damages', pol]] # get rid of geometry
    print(merged_df.dtypes)

    # merged_df = gpd.GeoDataFrame(merged_df) # reconvert to gdf after merge
    # merged_df[pol] = merged_df[pol] + merged_df[name]
    merged_df.fillna(0, inplace=True) # when adding nan to int, result was nan
    merged_df[pol] = merged_df.apply(lambda row: row[pol] + row['PM2.5 damages'], axis=1)
    # exit()
    # join with grid
    gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    merged_df = gdf.merge(merged_df,  how='inner', left_on=['i','j'], right_on = ['i','j'])[['i', 'j', pol, 'geometry']]
    # merged_df = merged_df[['i', 'j', pol, 'geometry']]
    # new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    # merged_df.to_file(new_name)
    return merged_df

""" normalized ld and hd """
def normalize_total_damages(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    # df = new_pmeds(fn, shpfn, pol, scc_name, scc)
    df = damages2(fn, shpfn, pol, name, points, scc_name, scc)

    # Get max total square
    df2 = damages3(fn, shpfn, pol, name, points, scc_name='Both', scc=[])
    # df2 = calculate_total2(fn, shpfn, pol, 'PM2.5T', 10, scc_name='Both', scc=[])
    max_v = max(df2[pol])

    df[pol] = df[pol]/max_v

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    # gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # # print(gdf.dtypes)
    #
    # # NOTE: Make sure dtypes of the columns in both tables are the same
    # # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # # x/j is 1 to 321, y/i is 1 to 291
    # # Keep only right or inner
    # merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[[pol, 'geometry']] # columns to keep, shape 599 by 2
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    df.to_file(new_name)

    # Add projection file
    from shutil import copy2
    copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

""" normalized facilities """
def normalize_total_damages2(fn, shpfn, pol, name, points, scc_name='', scc=[]):
    # df = new_pmeds(fn, shpfn, pol, scc_name, scc)
    df = damages(fn, points, shpfn, pol)

    # Get max total square
    df2 = damages3(fn, shpfn, pol, name, points, scc_name='Both', scc=[])
    # df2 = calculate_total2(fn, shpfn, pol, 'PM2.5T', 10, scc_name='Both', scc=[])
    max_v = max(df2[pol])
    name = 'PM2.5 damages'
    df[pol] = df[name]/max_v

    # i (int64) range 163 to 214
    # j (int64) range 115 to 149
    # df shape 599 by 20
    # print(df.shape)
    # print(df.head())
    # print(df['j'].min())
    # print(df['j'].max())
    print(df.shape)

    # print(df.sort_values(['pm'], ascending=[False]).head())
    # Join with grid shapefile
    # gdf = gpd.read_file(shpfn) # gdf shape 93411 by 4
    # # print(gdf.dtypes)
    #
    # # NOTE: Make sure dtypes of the columns in both tables are the same
    # # gdf['i'] = gdf['i'].map(lambda x: 291 - x) # switch which side i starts from
    # # gdf['j'] = gdf['j'].map(lambda x: x+1) # fix j being off by 1
    # # x/j is 1 to 321, y/i is 1 to 291
    # # Keep only right or inner
    # merged_df = gdf.merge(df,  how='right', left_on=['i','j'], right_on = ['i','j'])[[pol, 'geometry']] # columns to keep, shape 599 by 2
    # merged_df = gdf.merge(df, how='right', on=['i', 'j'])
    # print(merged_df.sort_values(['pm'], ascending=[False]).head())
    # print(merged_df.shape)
    new_name = save_name('MEDS/output/' + pol + fn.replace('MEDS/', '_') + scc_name, '.shp')
    df.to_file(new_name)

    # Add projection file
    # from shutil import copy2
    # copy2('MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))

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
