# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:49:06 2024

@author: fism
"""

import pandas as pd
import plotly.express as px
from dash import Dash, html

def heat_demand_plot(heat_35degC_demand, heat_65degC_demand):
    nHours = len(heat_35degC_demand)
    df = pd.DataFrame({
        'Time': range(nHours),
        'Heating Demand 35°C': heat_35degC_demand,
        'Heating Demand 65°C': heat_65degC_demand
    })
    df_melted = df.melt(id_vars='Time', value_vars=['Heating Demand 35°C', 'Heating Demand 65°C'],
                        var_name='Demand Type', value_name='Demand (kW)')
    fig = px.area(df_melted, x='Time', y='Demand (kW)', color='Demand Type',
                  title='Heat Demand 2022 NEST')
    return fig

def plot_power_generation(P_PV, P_imp, P_exp, df_input):
    df_plot = pd.DataFrame({
        'Time': range(len(df_input)),
        'PV Generation [kW]': [p / 1000 for p in P_PV],
        'Imported Power [kW]': [p / 1000 for p in P_imp],
        'Exported Power [kW]': [-p / 1000 for p in P_exp]
    })
    df_long = df_plot.melt(id_vars=['Time'], value_vars=['PV Generation [kW]', 'Imported Power [kW]', 'Exported Power [kW]'],
                           var_name='Type', value_name='Power [kW]')
    fig = px.line(df_long, x='Time', y='Power [kW]', color='Type',
                  title='PV Generation, Imported, and Exported Power Overview')
    return fig

def plot_component_sizes(S_PV, S_PV_max, S_ELY, S_ELY_max, S_C, S_C_max, S_FC, S_FC_max, S_TANK, S_TANK_max):
    components = ['S_PV', 'S_ELY', 'S_C', 'S_FC', 'S_TANK']
    sizes = [S_PV / 1000, S_ELY / 1000, S_C / 1000, S_FC / 1000, S_TANK / 3600000]
    max_sizes = [S_PV_max / 1000, S_ELY_max / 1000, S_C_max / 1000, S_FC_max / 1000, S_TANK_max / 3600000]
    df_components = pd.DataFrame({
        'Component': components,
        'Current Size': sizes,
        'Maximum Size': max_sizes
    })
    fig = px.bar(df_components, x='Component', y=['Current Size', 'Maximum Size'],
                 barmode='group', title='Component Sizes and Maximums')
    return fig

def plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input):
    nHours = len(df_input)
    df_HES = pd.DataFrame({
        'Time': range(nHours),
        'P_PV': [p / 1000 for p in P_PV],
        'P_ELY': [p / 1000 for p in P_ELY],
        'P_FC': [p / 1000 for p in P_FC],
        'E_TANK': [e / 3600000 for e in E_TANK]
    })
    fig = px.line(df_HES, x='Time', y=['P_PV', 'P_ELY', 'P_FC', 'E_TANK'],
                  labels={'value': 'Energy or Power [kW or MWh]', 'variable': 'Component'},
                  title='HESS Results Overview')
    return fig

def plot_costs_and_prices(cost_inst, cost_op, cost, cost_maint, df_input):
    costs = {'Installation Cost': cost_inst, 'Operation Cost': cost_op, 'Main Cost': cost_maint, 'Total Cost': cost}
    df_costs = pd.DataFrame(list(costs.items()), columns=['Cost Type', 'Amount'])
    df_prices = pd.DataFrame({
        'Time': range(len(df_input)),
        'Import Price': df_input['price_Eur_MWh'],
        'Export Price': df_input['Price_DayAhed']
    }).melt(id_vars='Time', var_name='Price Type', value_name='Price [€/MWh]')
    fig_costs = px.bar(df_costs, x='Cost Type', y='Amount', title='Costs Overview')
    fig_prices = px.line(df_prices, x='Time', y='Price [€/MWh]', color='Price Type', title='Electricity Prices Over Time')
    return fig_costs, fig_prices
