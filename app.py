# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import pandas as pd
import geopandas as gpd
from helpers import save_name, county_code, read_census_shape, read_facilities,\
agg_facilities_column, agg_facilities_allcols, create_grid, read_pmeds, agg_pmeds_column,\
combine_pmeds_facilities, facilities_to_pmeds, shape_pop, pop_emissions, shape_ethnicity,\
shape_income, shape_health, poverty_emissions

app = dash.Dash()

colors = {
    'background': '#FFFFFF',
    'text': '#7FDBFF'
}

df = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/' +
    '5d1ea79569ed194d432e56108a04d188/raw/' +
    'a9f9e8076b837d541398e999dcbac2b2826a81f8/'+
    'gdp-life-exp-2007.csv')
print(df.head())

gdf = gpd.read_file('out/combined/pov_emissionsFresno.shp')

markdown_text = '.'

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(
        children='CEC 2020',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    ),

    html.Div(children='Fresno Poverty Percentile to PM2.5 Emissions Scatter Plot', style={
        'textAlign': 'center',
        'color': colors['text']
    }),

    dcc.Graph(
        id='pov-emissions',
        figure={
            'data': [
                go.Scatter(
                    x=gdf[gdf['Tract_1'] == i]['Pov_pctl'],
                    y=gdf[gdf['Tract_1'] == i]['PM2.5T'],
                    # text=gdf[gdf['Tract_1'] == i]['Tract_1'],
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                    name=i
                ) for i in gdf.Tract_1.unique()
            ],
            'layout': go.Layout(
                xaxis={'type': 'log', 'title': 'Poverty Percentile (Statewide)'},
                yaxis={'title': 'PM2.5 (Raw)'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                # title='\nFresno Poverty Percentile to PM2.5 Emissions Scatter Plot'
            )
        }
    ),

    dcc.Graph(
        id='pov-emissions-map',
        figure={
            'data': [
                go.Choropleth(
                    locationmode = "ISO-3",
                    locations = gdf[gdf['Tract_1'] == i]['Tract_1'],
                    z = gdf[gdf['Tract_1'] == i]['PM2.5T'],
                    text = gdf[gdf['Tract_1'] == i]['Tract_1'],
                    colorscale = [[0,'rgb(0, 0, 0)'],[1,'rgb(0, 0, 0)']],
                    autocolorscale = False,
                    showscale = False,
                    # geo = 'geo2'
                    geo = gdf[gdf['Tract_1'] == i]['geometry'].to_json()
                ) for i in gdf.Tract_1.unique()
            ]
            # 'layout': go.Layout(
            #     xaxis={'type': 'log', 'title': 'Poverty Percentile (Statewide)'},
            #     yaxis={'title': 'PM2.5 (Raw)'},
            #     margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            #     legend={'x': 0, 'y': 1},
            #     hovermode='closest',
            #     title=''
            # )
        }
    ),

    dcc.Graph(
        id='example-graph-2',
        figure={
            'data': [
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'bar', 'name': 'SF'},
                {'x': [1, 2, 3], 'y': [2, 4, 5], 'type': 'bar', 'name': u'Montréal'},
            ],
            'layout': {
                'plot_bgcolor': colors['background'],
                'paper_bgcolor': colors['background'],
                'font': {
                    'color': colors['text']
                }
            }
        }
    ),

    dcc.Graph(
        id='life-exp-vs-gdp',
        figure={
            'data': [
                go.Scatter(
                    x=df[df['continent'] == i]['gdp per capita'],
                    y=df[df['continent'] == i]['life expectancy'],
                    text=df[df['continent'] == i]['country'],
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size': 15,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                    name=i
                ) for i in df.continent.unique()
            ],
            'layout': go.Layout(
                xaxis={'type': 'log', 'title': 'GDP Per Capita'},
                yaxis={'title': 'Life Expectancy'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    ),

    dcc.Markdown(children=markdown_text)
])

if __name__ == '__main__':
    app.run_server(debug=True)
