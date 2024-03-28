# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:32:04 2024

Dashboard to visualize the results from pv_sizingprob_WHR.

You need to install plotly and dash before running the code:

1. plotly may be installed using pip: https://plotly.com/python/getting-started/
    $ pip install plotly==5.20.0
    or conda:
    $ conda install -c plotly plotly=5.20.0

This package contains everything you need to write figures to standalone HTML files.

2. In your terminal, install dash: https://dash.plotly.com/installation
    $ pip install dash
    or conda:
    $ conda install conda-forge::dash

Your default web browser should automatically open to http://127.0.0.1:8050/, 
displaying the dashboard. 
If it doesn't open automatically, you can manually open your web browser and 
enter the address http://127.0.0.1:8050/ into the address bar to view 
the dashboard.

@author: fism
"""

# Import necessary components from Dash and Plotly
from dash import Dash, dcc, html, Input, Output, dash_table
import os
import threading
import webbrowser

# Path to pv_sizingprob_WHR
path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\PSWHR\PSWHR\PSWHR'
os.chdir(path)
# from pv_sizingprob_WHR import *
from pv_sizingprob_WHR import plot_power_generation, P_PV, P_imp, P_exp, df_input, nHours
from pv_sizingprob_WHR import costs_pie_chart, all_costs
from pv_sizingprob_WHR import plot_WHR, results
                                    


app = Dash(__name__, suppress_callback_exceptions=True)
app.title = "PSWHR Dashboard"

# Define the app layout with navigation links
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div([
        dcc.Link('Demand | ', href='/demand'),
        dcc.Link('Sizes | ', href='/sizes'),
        dcc.Link('Operation | ', href='/operation'),
        dcc.Link('WHR | ', href='/WHR'),
        dcc.Link('Costs', href='/costs'),
    ], className='row'),
    html.Div(id='page-content')
])

# Layout for the Demand page
demand_layout = html.Div([
    html.H1('Power Demand and inputs'),
    html.P('The electrical power demand comes from the Tesla Supecharger Station in Rubigen, Switzerland.'),
    # Add a Graph component where your Plotly graph will be rendered
    dcc.Graph(id='power-generation-graph')
])

# Layout for the Sizes page
sizes_layout = html.Div([
    html.H1('Hess and Battery Sizes'),
    html.P('Information related to sizes.')
])

# Layout for the Operation page
operation_layout = html.Div([
    html.H1('Operation Page'),
    html.P('Information related to operation.')
])

# Layout for the Sizes page
WHR_layout = html.Div([
    html.H1('Waste Heat Recovery'),
    html.P('Heat is recovered from both the Electrolyzer and the Fuel Cell.'),
    dcc.Graph(id='WHR', figure=plot_WHR(results))
])

costs_layout = html.Div([
    html.H1('System Costs'),
    html.Div([
        dcc.Graph(
            id='costs-pie-chart', figure=costs_pie_chart(all_costs), 
            style={'width': '50%', 'display': 'inline-block'}),
        dash_table.DataTable(
            id='costs-table',
            style_table={'width': '50%', 'display': 'inline-block', 'paddingLeft': '20px'}
        ),
    ], style={'display': 'flex'})
])

# Callback to switch pages and update the page content
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/demand':
        return demand_layout
    elif pathname == '/sizes':
        return sizes_layout
    elif pathname == '/operation':
        return operation_layout
    elif pathname == '/WHR':
        return WHR_layout
    elif pathname == '/costs':
        return  costs_layout
    else:
        return '404'

@app.callback(Output('power-generation-graph', 'figure'), [Input('url', 'pathname')])
def update_graph(pathname):
    if pathname == '/demand':
        return plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours)
    return {}

@app.callback(
    [Output('costs-table', 'data'), Output('costs-table', 'columns')],
    [Input('url', 'pathname')]
)
def update_table(pathname):
    if pathname == '/costs':
        data = [{"Cost Component": key, "Value in kCHF": value} for key, value in all_costs.items()]
        columns = [{"name": i, "id": i} for i in data[0].keys()]
        return data, columns
    return [], []

# Function to open web browser at the app URL
def open_browser():
      webbrowser.open_new("http://127.0.0.1:8050/")

# Function to run the server
def run_dash_app():
    app.run_server(debug=False)  # Set debug=False to prevent opening two tabs

# Run the app in a separate thread
if __name__ == '__main__':
    threading.Thread(target=run_dash_app).start()
    # Wait for the server to start up and then open the browser
    threading.Timer(1.5, open_browser).start()
