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
        area: 'blockgroup', 'block', 'pop'
        county: 'LA', 'Fresno'
    """
    if area == 'blockgroup':
        area = gpd.read_file('data/ca_blockgroup/tl_2016_06_bg.shp')
    elif area == 'block':
        area = gpd.read_file('data/ca_block/tl_2016_06_tabblock10.shp')
    elif area == 'pop': # block level
        area = gpd.read_file('data/fresno_pop/fresnopop.shp')
    else:
        raise ValueError('check read_census_shape for valid areas')

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

def agg_facilities_column(area, county, column, output='percent', damages=False, save=False):
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
        df = df[df['SCC'].isin([302,303,304,305,306,307,309,310,311,312,313,314,315,316,318,319,321,408,413,420,414,508,513,514,520,613,614,617,620,717])]

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
        print("Saved " + new_name)

    return merged_df

def combine_pmeds_facilities(county, column, scc='', output='percent', updated=True, damages=False, save=False):
    """
    Return geodataframe of column added from PMEDS and facilities for each grid cell
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'TOGT'
        scc: 'LD', 'HD', '' for all
        output: 'percent' for percent of total that each grid cell is, 'actual'
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
        print("Saved " + new_name)

    return df

def facilities_to_pmeds(county, column, updated=True, save=False):
    """
    Return geodataframe of percent that HD, LD, all vehicles (Both), and I are of total of column for each grid cell
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'TOGT'
        updated: True to update to 2015 values
        save: Whether to save frame to shapefile
    """
    # Calculate raw PMEDS frames
    pmeds_ld = agg_pmeds_column(county, column, 'LD', 'actual', updated)
    pmeds_hd = agg_pmeds_column(county, column, 'HD', 'actual', updated)
    pmeds = agg_pmeds_column(county, column, '', 'actual', updated)
    print(pmeds.shape)
    # Calculate raw facilities frame
    facilities = agg_facilities_column('grid', county, column, 'actual')
    print(facilities.shape)
    total = pd.concat([pmeds, facilities])
    # Add up individual squares' pm
    total[column] = total.groupby(['i', 'j'])[column].transform('sum') # Add each row containing same i, j, but keep dupes
    total = total.drop_duplicates(['i', 'j']) # Get rid of dupes
    print(total.shape)
    total = total[total[column] != 0]
    print(total.shape)

    # Industry
    total = total.merge(facilities,  how='left', on=['i','j'])[[column + '_x', column + '_y', 'geometry_x', 'i', 'j']] # columns to keep
    total['I'] = total.apply(lambda row: (row[column+'_y']/row[column+'_x'])*100, axis=1) # facilities/total
    total.rename(columns = {'geometry_x':'geometry', column + '_x': column}, inplace = True)
    total = total.drop(column + '_y', 1)
    total = gpd.GeoDataFrame(total)

    # LD
    total = total.merge(pmeds_ld,  how='left', on=['i','j'])[[column + '_x', column + '_y', 'geometry_x', 'i', 'j', 'I']] # columns to keep
    total['LD'] = total.apply(lambda row: (row[column+'_y']/row[column+'_x'])*100, axis=1) # ld/total
    total.rename(columns = {'geometry_x':'geometry', column + '_x': column}, inplace = True)
    total = total.drop(column + '_y', 1)
    total = gpd.GeoDataFrame(total)

    # HD
    total = total.merge(pmeds_hd,  how='left', on=['i','j'])[[column + '_x', column + '_y', 'geometry_x', 'i', 'j', 'I', 'LD']] # columns to keep
    total['HD'] = total.apply(lambda row: (row[column+'_y']/row[column+'_x'])*100, axis=1) # hd/total
    total.rename(columns = {'geometry_x':'geometry', column + '_x': column}, inplace = True)
    total = total.drop(column + '_y', 1)
    total = gpd.GeoDataFrame(total)

    # Both LD and HD
    total = total.merge(pmeds,  how='left', on=['i','j'])[[column + '_x', column + '_y', 'geometry_x', 'i', 'j', 'I', 'LD', 'HD']] # columns to keep
    total['Both'] = total.apply(lambda row: (row[column+'_y']/row[column+'_x'])*100, axis=1) # both/total
    total.rename(columns = {'geometry_x':'geometry', column + '_x': column}, inplace = True)
    total = total.drop(column + '_y', 1)
    total = gpd.GeoDataFrame(total)

    def dominant_percentage(row):
        max_column = max(row['I'], row['LD'], row['HD'])
        if max_column == row['I']:
            return 'I'
        if max_column == row['LD']:
            return 'LD'
        return 'HD'

    total['Dominant'] = total.apply(dominant_percentage, axis=1)

    print(total.head())
    print(total.dtypes)

    if save:
        new_name = save_name('out/combined/' + county + column, '.shp')
        total.to_file(new_name)
        # Add projection file
        from shutil import copy2
        copy2('data/MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))
        print("Saved " + new_name)

    return total

def shape_pop(area, save=False):
    """
    Return geodataframe of population in each polygon of area added up. Uses centroids of original population blocks to match to new area.
        area: 'blockgroup', 'grid'
        county: 'LA', 'Fresno'
        save: Whether to save frame to shapefile
    """
    if area == 'blockgroup':
        pop = read_census_shape('pop')
        pop['geometry'] = pop['geometry'].centroid # set all blocks to points
        gdf = read_census_shape('blockgroup')
        total = gpd.sjoin(pop, gdf, op='within')
        total = total.groupby('index_right')['POP10'].sum().to_frame() # need to convert Series to DF
        total = gdf.join(total).dropna()[['GEOID', 'POP10', 'geometry']]
    else:
        pop = read_census_shape('pop').to_crs(epsg=4326)
        pop['geometry'] = pop['geometry'].centroid # set all blocks to points
        gdf = gpd.read_file('data/MEDS/grid.shp').to_crs(epsg=4326) # gdf shape 93411 by 4
        # Gives table of pop with column 'index_right' being which polygon each falls in
        total = gpd.sjoin(pop, gdf, op='within')
        total = total.groupby('index_right')['POP10'].sum().to_frame() # need to convert Series to DF
        total = gdf.join(total).dropna()[['i', 'j', 'POP10', 'geometry']]

    if save:
        new_name = save_name('out/combined/pop_' + area, '.shp')
        total.to_file(new_name)
        # No projection needed since coordinate system defined
        print("Saved " + new_name)

    # TODO: pop['density'] = pop['POP10']/pop['geometry'].area

    return total

def pop_emissions(area, county, column, save=False):
    """
    Return geodataframe of population in each polygon of area added up and multiplied by pollution
        area: 'blockgroup', 'grid'
        county: 'LA', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'TOGT'
        save: Whether to save frame to shapefile
    """
    if area == 'blockgroup':
        pop = shape_pop('blockgroup').to_crs(epsg=4326) # [['GEOID', 'POP10']] # get rid of geometry for join
        em = combine_pmeds_facilities(county, column, output='actual', updated=True).to_crs(epsg=4326)
        em['geometry'] = em['geometry'].centroid # set all grid cells to points
        total = gpd.sjoin(em, pop, op='within')
        total = total.groupby('index_right')[column].sum().to_frame() # need to convert Series to DF
        total = pop.join(total).dropna() # [['i', 'j', 'POP10', 'geometry']]
        total['POP10'] = total[column] * total['POP10']
        total = total[['GEOID', 'POP10', 'geometry']]
    else:
        pop = shape_pop('grid')[['i', 'j', 'POP10']] # get rid of geometry for join
        em = combine_pmeds_facilities(county, column, output='actual', updated=True)

        total = em.merge(pop,  how='outer', on=['i', 'j'])[['i', 'j', column, 'POP10', 'geometry']] # can't hash on geometry
        total['POP10'] = total[column] * total['POP10']
        total = total[['i', 'j', 'POP10', 'geometry']]

    if save:
        new_name = save_name('out/combined/' + county + column, '.shp')
        total.to_file(new_name)
        if area == 'grid':
            # Add projection file
            from shutil import copy2
            copy2('data/MEDS/_mm5_sphere_.prj', new_name.replace('.shp', '.prj'))
        print("Saved " + new_name)

    return total

def shape_ethnicity(area, county, save=False):
    """
    Return geodataframe of percent that each ethnicity is in shape
        area: 'block', 'grid'
        county: 'Los Angeles', 'Fresno'
        save: Whether to save frame to shapefile
    """
    gdf = gpd.read_file('data/ces3shp/CES3Results.shp').to_crs(epsg=4326)
    gdf = gdf[gdf['County'] == county][['Hispanic__', 'White____', 'African_Am', 'Native_Ame', 'Asian_Amer', 'Other____', 'geometry']]

    # Convert fields from object to float
    gdf = gdf[gdf['Hispanic__'] != '<Null>']
    gdf['Hispanic__'] = gdf['Hispanic__'].astype(str).astype(float)
    gdf = gdf[gdf['White____'] != '<Null>']
    gdf['White____'] = gdf['White____'].astype(str).astype(float)
    gdf = gdf[gdf['African_Am'] != '<Null>']
    gdf['African_Am'] = gdf['African_Am'].astype(str).astype(float)
    gdf = gdf[gdf['Native_Ame'] != '<Null>']
    gdf['Native_Ame'] = gdf['Native_Ame'].astype(str).astype(float)
    gdf = gdf[gdf['Asian_Amer'] != '<Null>']
    gdf['Asian_Amer'] = gdf['Asian_Amer'].astype(str).astype(float)
    gdf = gdf[gdf['Other____'] != '<Null>']
    gdf['Other____'] = gdf['Other____'].astype(str).astype(float)

    if area == 'grid':
        # Squares get value of whatever block their centroid is in
        grid = gpd.read_file('data/MEDS/grid.shp').to_crs(epsg=4326) # gdf shape 93411 by 4
        grid['geometry'] = grid['geometry'].centroid # set all blocks to points
        # Gives table of pop with column 'index_right' being which polygon each falls in
        gdf = gpd.sjoin(grid, gdf, op='within')[['Hispanic__', 'White____', 'African_Am', 'Native_Ame', 'Asian_Amer', 'Other____', 'i', 'j']]
        grid = gpd.read_file('data/MEDS/grid.shp').to_crs(epsg=4326) # gdf shape 93411 by 4
        gdf = grid.merge(gdf,  how='right', on=['i','j'])[['Hispanic__', 'White____', 'African_Am', 'Native_Ame', 'Asian_Amer', 'Other____', 'i', 'j', 'geometry']]

    def dominant_percentage(row):
        max_column = max(row['Hispanic__'], row['White____'], row['African_Am'], row['Native_Ame'], row['Asian_Amer'], row['Other____'])
        if max_column == row['Hispanic__']:
            return 'Hispanic'
        if max_column == row['White____']:
            return 'White'
        if max_column == row['African_Am']:
            return 'African_Am'
        if max_column == row['Native_Ame']:
            return 'Native_Am'
        if max_column == row['Asian_Amer']:
            return 'Asian_Am'
        return 'Other'

    gdf['Dominant'] = gdf.apply(dominant_percentage, axis=1)

    print(gdf.head())

    if save:
        new_name = save_name('out/combined/eth' + county, '.shp')
        gdf.to_file(new_name)
        print("Saved " + new_name)

    return gdf

def shape_income(county, save=False):
    """
    Return geodataframe of percentile of poverty in shape
        area: 'block', 'grid'
        county: 'Los Angeles', 'Fresno'
        save: Whether to save frame to shapefile
    """
    gdf = gpd.read_file('data/ces3shp/CES3Results.shp').to_crs(epsg=4326)
    gdf = gdf[gdf['County'] == county][['Pov_pctl', 'geometry', 'Tract_1']]

    if save:
        new_name = save_name('out/combined/inc' + county, '.shp')
        gdf.to_file(new_name)
        print("Saved " + new_name)

    return gdf

def shape_health(county, save=False):
    """
    Return geodataframe of percentile of health effects in shape
        area: 'block', 'grid'
        county: 'Los Angeles', 'Fresno'
        save: Whether to save frame to shapefile
    """
    gdf = gpd.read_file('data/ces3shp/CES3Results.shp').to_crs(epsg=4326)
    gdf = gdf[gdf['County'] == county][['CVD_pctl', 'LBW_pctl', 'Asthma_Pct', 'geometry']]

    if save:
        new_name = save_name('out/combined/health' + county, '.shp')
        gdf.to_file(new_name)
        print("Saved " + new_name)

    return gdf

def poverty_emissions(county, column, save=False):
    """
        county: 'Los Angeles', 'Fresno'
        column: 'PM2.5T', 'PMT', 'SOXT', 'NOXT', 'COT', 'TOGT'
        save: Whether to save frame to shapefile
    """
    pov = shape_income(county)
    # county = 'LA' if county == 'Los Angeles'
    em = combine_pmeds_facilities(county, column, output='actual', updated=True).to_crs(epsg=4326)
    em['geometry'] = em['geometry'].centroid # set all grid cells to points
    total = gpd.sjoin(em, pov, op='within')
    total = total.groupby('index_right')[column].sum().to_frame() # need to convert Series to DF
    total = pov.join(total).dropna()

    if save:
        new_name = save_name('out/combined/pov_emissions' + county, '.shp')
        total.to_file(new_name)
        print("Saved " + new_name)

    return total
