# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:25:19 2024

@author: fism
"""

import matplotlib.pyplot as plt
import seaborn as sns

# Set Seaborn style
sns.set(style="whitegrid")
palette = sns.color_palette("husl", 8)  # You can choose a different palette if needed

def plot_power_generation(P_PV, P_imp, P_exp, df_input):
    # Plot power generation, imported, and exported power
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))  # Use a single axis
    ax.plot(range(len(df_input)), [P_PV[t]  / 1000 for t in range(len(df_input))], label='PV Generation', color='orange')
    ax.plot(range(len(df_input)), [P_imp[t] / 1000 for t in range(len(df_input))], label='Imported Power', color='blue')
    # Plot P_exp as negative values
    ax.plot(range(len(df_input)), [-P_exp[t] / 1000 for t in range(len(df_input))], label='Exported Power', color='green')
    
    ax.set_ylabel('Power [kW]')
    ax.legend(loc='upper left')
    ax.set_title('PV Generation, Imported, and Exported Power')
    ax.set_xlabel('Time [h]')  # Set the x-label here

    plt.show()

    
def plot_component_sizes(Area_PV, Area_PV_max, S_ELY, S_FC, S_TANK, S_ELY_max, S_FC_max, S_TANK_max):
    # Plot sizes of components (ELY, FC, TANK)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
    
    components_axis1 = ['Area_PV', 'Area_PV_max']
    component_sizes_axis1 = [Area_PV, Area_PV_max]
    colors_axis1 = ['mediumturquoise', 'aquamarine']

    ax1.bar(components_axis1, component_sizes_axis1, color=colors_axis1, width=0.7)
    #ax1.set_xlabel('Components')
    ax1.set_ylabel('Size (Squaremeters)')
    
    components_axis2 = ['S_ELY', 'S_ELY_max', 'S_FC', 'S_FC_max']
    component_sizes_axis2 = [S_ELY / 1000, S_ELY_max / 1000, S_FC / 1000, S_FC_max / 1000]
    colors_axis2 = ['blue', 'lightblue', 'green', 'lightgreen']
    
    ax2.bar(components_axis2, component_sizes_axis2, color=colors_axis2, width=0.7)
    #ax2.set_xlabel('Components')
    ax2.set_ylabel('Size (Kilowatts)')

    components_axis3 = ['S_TANK', 'S_TANK_max']
    component_sizes_axis3 = [S_TANK / 3600000, S_TANK_max / 3600000]
    colors_axis3 = ['orange', 'lightcoral']

    ax3.bar(components_axis3, component_sizes_axis3, color=colors_axis3, width=0.7)
    #ax3.set_xlabel('Components')
    ax3.set_ylabel('Size (Kilowatt-hours)')

    plt.tight_layout()
    plt.show()


def plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input):
    """
    Plot main results including the size and operation of the electrolyzer (ELY),
    fuel cell (FC), and energy and TANK size over time.

    Parameters:
    - P_PV: list, hourly PV power generation values
    - P_ELY: list, hourly electrolyzer power values
    - S_ELY: float, electrolyzer size
    - P_FC: list, hourly fuel cell power values
    - S_FC: float, fuel cell size
    - E_TANK: list, hourly tank energy values
    - S_TANK: float, tank size
    - S_TANK_max: float, maximal tank size
    - df_input: DataFrame, input data
    """
    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(20, 10))

    # Plot 1: Size and operation of the ELY
    ax1.plot(range(len(df_input)), [P_ELY[t] / 1000 for t in range(len(df_input))], label='ELY Electrolyzer (kW)', color=palette[1])
    ax1.plot(range(len(df_input)), [S_ELY / 1000] * len(df_input), label='ELY size (kW)', linestyle='--', color=palette[3])
    ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Power [kW]')
    ax1.set_ylim([0, S_ELY_max * 1.1 / 1000])
    ax1.legend(loc='upper right')

    # Plot 2: Size and operation of the FC
    ax2.plot(range(len(df_input)), [P_FC[t] / 1000 for t in range(len(df_input))], label='Fuel Cell (kW)', color=palette[0])
    ax2.plot(range(len(df_input)), [S_FC / 1000] * len(df_input), label='FC size (kW)', linestyle=':', color=palette[2])
    ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Power [kW]')
    ax2.set_ylim([0, S_FC_max * 1.1 / 1000])
    ax2.legend(loc='upper right')

    # Plot 3: Energy and TANK Size
    ax3.plot(range(len(df_input)), [E_TANK[t] / 3600000 for t in range(len(df_input))], label='TANK (kWh)', color=palette[5])
    ax3.plot(range(len(df_input)), [S_TANK / 3600000] * len(df_input), label='TANK size (kWh)', linestyle='-.', color=palette[6])
    ax3.set_xlabel('Time [h]')
    ax3.set_ylabel('Energy [kWh]')
    ax3.set_ylim([0, S_TANK_max * 1.1 / 1000 / 3600])
    ax3.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

def plot_costs_and_prices(cost_inst, cost_op, cost, cost_maint, cost_startup, df_input):
    """
    Plot installation cost, operation cost, total cost, and electricity prices.

    Parameters:
    - cost_inst: LinExpr object of gurobipy module, installation cost
    - cost_op: LinExpr object of gurobipy module, operation cost
    - cost: LinExpr object of gurobipy module, total cost
    - cost_main: LinExpr object of gurobipy module, main cost
    - cost_startup: LinExpr object of gurobipy module, startup cost
    - df_input: DataFrame, input data
    """

    # Set Seaborn style
    sns.set(style="whitegrid")
    palette = sns.color_palette("pastel")

    # Plotting
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 10))

    # Plot electricity prices
    ax1.plot(range(len(df_input)), df_input['price_Eur_MWh'], label='Import Price (EUR/MWh)', color=palette[0])
    ax1.plot(range(len(df_input)), df_input['Price_DayAhed'], label='Export Price (EUR/MWh)', color=palette[1])
    ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Price [EUR/Wh]')
    ax1.legend()

    # Plot installation cost, operation cost, main cost, startup cost, and total cost
    bar_width = 0.4
    bar1 = ax2.bar(0, cost_inst / 1000, bar_width, label='Installation Cost (k€/y)', color=palette[2])
    bar2 = ax2.bar(1, cost_op   / 1000, bar_width, label='Operation Cost (k€/y)', color=palette[3])
    bar3 = ax2.bar(2, cost_maint / 1000, bar_width, label='Main Cost (k€/y)', color=palette[4])
    bar4 = ax2.bar(3, cost_startup / 1000, bar_width, label='Startup Cost (k€/y)', color=palette[5])
    bar5 = ax2.bar(4, cost      / 1000, bar_width, label='Total Cost (k€/y)', color=palette[6])

    ax2.set_ylabel('Cost [k€/y]')
    ax2.legend()
    ax2.set_xticks([0, 1, 2, 3, 4])
    ax2.set_xticklabels(['Installation Cost', 'Operation Cost', 'Maintenance Cost', 'Startup Cost', 'Total Cost'])

    plt.tight_layout()
    plt.show()

