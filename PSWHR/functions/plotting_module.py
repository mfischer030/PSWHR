# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:25:19 2024

@author: fism
"""
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set Seaborn style
sns.set(style="whitegrid")
palette = sns.color_palette("husl", 8)  # You can choose a different palette if needed

# Increase font sizes for better readability on PowerPoint slides
plt.rcParams.update({'font.size': 24})  # Adjust global font size
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['legend.fontsize'] = 24

def heat_demand_plot(heat_35degC_demand,heat_65degC_demand):
    # Creating figure and subplots for heat demand
    fig, axs = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    fig.suptitle('Heat Demand 2022 NEST', fontsize=28)

    # Plot Heating demand at 35°C
    axs[0].fill_between(range(len(heat_35degC_demand)), heat_35degC_demand, label='Heating Demand 35°C', color='orange', alpha=0.4)
    axs[0].set_ylabel('Demand (kW)', fontsize=24)
    axs[0].legend(loc='upper right', fontsize=24)

    # Plot DHW demand at 65°C
    axs[1].fill_between(range(len(heat_65degC_demand)), heat_65degC_demand, label='Heat Demand 65°C', color='red', alpha=0.4)
    axs[1].set_ylabel('Demand (kW)', fontsize=24)
    axs[1].legend(loc='upper right', fontsize=24)

    # Setting the xlabel for the last subplot
    axs[1].set_xlabel('Time in Hours', fontsize=24)

    plt.tight_layout()
    plt.show()

# def plot_power_generation(P_PV, P_imp, P_exp, df_input):
#     # Plot power generation, imported, and exported power
#     fig, ax = plt.subplots(1, 1, figsize=(20, 10))  # Use a single axis
#     ax.plot(range(len(df_input)), [P_PV[t]  / 1000 for t in range(len(df_input))], label='PV Generation', color='orange')
#     ax.plot(range(len(df_input)), [P_imp[t] / 1000 for t in range(len(df_input))], label='Imported Power', color='blue')
#     # Plot P_exp as negative values
#     ax.plot(range(len(df_input)), [-P_exp[t] / 1000 for t in range(len(df_input))], label='Exported Power', color='green')
    
#     ax.set_ylabel('Power [kW]')
#     ax.legend(loc='upper left')
#     ax.set_title('PV Generation, Imported, and Exported Power')
#     ax.set_xlabel('Time [h]')  # Set the x-label here

#     plt.show()

def plot_power_generation(P_PV, P_imp, P_exp, df_input):
    # Assuming df_input has an index that can serve as the x-axis (time)
    x_values = df_input.index  # or range(len(df_input)) if there's no specific index
    
    # Create a DataFrame for Plotly Express
    df_plot = pd.DataFrame({
        'Time': x_values,
        'PV Generation [kW]': [p / 1000 for p in P_PV],
        'Imported Power [kW]': [p / 1000 for p in P_imp],
        'Exported Power [kW]': [-p / 1000 for p in P_exp]
    })

    # Melt the DataFrame to long format for easier plotting with Plotly Express
    df_long = pd.melt(df_plot, id_vars=['Time'], value_vars=['PV Generation [kW]', 'Imported Power [kW]', 'Exported Power [kW]'],
                      var_name='Type', value_name='Power')

    # Create the plot
    fig = px.line(df_long, x='Time', y='Power', color='Type',
                  labels={'Power': 'Power [kW]', 'Time': 'Time [h]'},
                  color_discrete_map={
                      'PV Generation [kW]': 'orange',
                      'Imported Power [kW]': 'magenta',
                      'Exported Power [kW]': 'green'
                  })

    # Update layout for aesthetics
    fig.update_layout(title='PV Generation, Imported, and Exported Power Overview',
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    
    fig.write_html('first_figure.html', auto_open=True)
    
# def plot_power_generation(P_PV, P_imp, P_exp, df_input):
#     # Plot power generation, imported, and exported power
#     fig, ax = plt.subplots(1, 1, figsize=(20, 6))  # Use a single axis with a larger figsize for readability
    
#     # Plot and fill PV Generation
#     ax.fill_between(range(len(df_input)), [P_PV[t] / 1000 for t in range(len(df_input))], color='orange', alpha=0.3)
#     ax.plot(range(len(df_input)), [P_PV[t] / 1000 for t in range(len(df_input))], label='PV Generation', color='orange')
    
#     # Plot and fill Imported Power
#     ax.fill_between(range(len(df_input)), [P_imp[t] / 1000 for t in range(len(df_input))], color='magenta', alpha=0.3)
#     ax.plot(range(len(df_input)), [P_imp[t] / 1000 for t in range(len(df_input))], label='Imported Power', color='magenta')
    
#     # Plot and fill Exported Power as negative values
#     ax.fill_between(range(len(df_input)), [-P_exp[t] / 1000 for t in range(len(df_input))], color='green', alpha=0.3)
#     ax.plot(range(len(df_input)), [-P_exp[t] / 1000 for t in range(len(df_input))], label='Exported Power', color='green')
    
#     ax.set_ylabel('Power [kW]')  # Increase the font size for the y-axis label
#     ax.set_xlabel('Time [h]')  # Increase the font size for the x-axis label
#     ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17),
#           fancybox=False, shadow=False, ncol=5)
#     #ax.legend(loc='upper right')  # Increase the font size for the legend
#     #ax.set_title('PV Generation, Imported, and Exported Power')  
    
#     # Increase the font size for the tick labels
#     ax.tick_params(axis='both', which='major')
    
#     plt.show()


# def plot_component_sizes(S_PV, S_PV_max, S_ELY, S_FC, S_TANK, S_ELY_max, S_FC_max, S_TANK_max):
#     # Plot sizes of components (PV, ELY, FC, TANK)
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 8))
    
#     # Integrating S_PV and S_PV_max into the first subplot (now ax1)
#     components_axis1 = ['S_PV', 'S_PV_max', 'S_ELY', 'S_ELY_max', 'S_FC', 'S_FC_max']
#     component_sizes_axis1 = [S_PV /1000, S_PV_max / 1000, S_ELY / 1000, S_ELY_max / 1000, S_FC / 1000, S_FC_max / 1000]
#     colors_axis1 = ['mediumturquoise', 'aquamarine', 'blue', 'lightblue', 'green', 'lightgreen']
    
#     ax1.bar(components_axis1, component_sizes_axis1, color=colors_axis1, width=0.7)
#     ax1.set_ylabel('Size (Kilowatts)')
    
#     # Second subplot (now ax2) remains unchanged, for TANK sizes
#     components_axis2 = ['S_TANK', 'S_TANK_max']
#     component_sizes_axis2 = [S_TANK / 3600000, S_TANK_max / 3600000]  # Converting to Kilowatt-hours
#     colors_axis2 = ['orange', 'lightcoral']
    
#     ax2.bar(components_axis2, component_sizes_axis2, color=colors_axis2, width=0.7)
#     ax2.set_ylabel('Size (Kilowatt-hours)')

#     plt.tight_layout()
#     plt.show()

def plot_component_sizes(S_PV, S_PV_max, S_ELY, S_ELY_max, S_C, S_C_max, S_FC, S_FC_max, S_TANK, S_TANK_max):
    # Setting the Seaborn color palette
    sns.set_theme(style="whitegrid")
    colors = sns.color_palette("deep", 5)  # Using the 'deep' palette for 5 different colors

    # First plot for S_PV, S_ELY, S_C, S_FC
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 8))
    components_kw = ['S_PV', 'S_ELY', 'S_C', 'S_FC']
    component_sizes_kw = [S_PV / 1000, S_ELY / 1000, S_C /1000, S_FC / 1000]  # Converting to kW
    max_values_kw = [S_PV_max / 1000, S_ELY_max / 1000, S_C_max / 1000, S_FC_max / 1000]
    
    ax0.bar(components_kw, component_sizes_kw, color=colors[:4], width=0.4)
    for i, max_val in enumerate(max_values_kw):
        ax0.axhline(y=max_val, color=colors[i], linestyle='--', label=f'{components_kw[i]} Max')

    ax0.set_ylabel('Size (Kilowatts)', fontsize=24)
    ax0.tick_params(axis='both', which='major', labelsize=24)
    
    # Second plot for S_TANK
    components_MWh = ['S_TANK']
    component_sizes_MWh = [S_TANK / (3600000*1000)]  # Converting to MWh
    color_MWh = colors[4]  # Using the fifth color for S_TANK

    ax1.bar(components_MWh, component_sizes_MWh, color=color_MWh, width=0.4)
    ax1.axhline(y=S_TANK_max / (3600000*1000), color=color_MWh, linestyle='--', label='S_TANK Max')

    ax1.set_ylabel('Size (Megawatt-hours)', fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)

    plt.tight_layout()
    plt.show()
    
# def plot_component_sizes(Area_PV, Area_PV_max, S_ELY, S_FC, S_TANK, S_ELY_max, S_FC_max, S_TANK_max):
#     # Plot sizes of components (ELY, FC, TANK)
#     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
    
#     components_axis1 = ['Area_PV', 'Area_PV_max']
#     component_sizes_axis1 = [Area_PV, Area_PV_max]
#     colors_axis1 = ['mediumturquoise', 'aquamarine']

#     ax1.bar(components_axis1, component_sizes_axis1, color=colors_axis1, width=0.7)
#     #ax1.set_xlabel('Components')
#     ax1.set_ylabel('Size (Squaremeters)')
    
#     components_axis2 = ['S_ELY', 'S_ELY_max', 'S_FC', 'S_FC_max']
#     component_sizes_axis2 = [S_ELY / 1000, S_ELY_max / 1000, S_FC / 1000, S_FC_max / 1000]
#     colors_axis2 = ['blue', 'lightblue', 'green', 'lightgreen']
    
#     ax2.bar(components_axis2, component_sizes_axis2, color=colors_axis2, width=0.7)
#     #ax2.set_xlabel('Components')
#     ax2.set_ylabel('Size (Kilowatts)')

#     components_axis3 = ['S_TANK', 'S_TANK_max']
#     component_sizes_axis3 = [S_TANK / 3600000, S_TANK_max / 3600000]
#     colors_axis3 = ['orange', 'lightcoral']

#     ax3.bar(components_axis3, component_sizes_axis3, color=colors_axis3, width=0.7)
#     #ax3.set_xlabel('Components')
#     ax3.set_ylabel('Size (Kilowatt-hours)')

#     plt.tight_layout()
#     plt.show()


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
    ax1.plot(range(len(df_input)), [P_ELY[t] / 1000 for t in range(len(df_input))], label='ELY', color=palette[1])
    ax1.plot(range(len(df_input)), [S_ELY / 1000] * len(df_input), label='ELY size', linestyle='--', color=palette[3])
    # ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Power [kW]')
    # ax1.set_ylim([0, S_ELY_max * 1.1 / 1000])
    ax1.set_ylim([0, S_ELY / 1000 * 1.1])
    ax1.legend(loc='upper right')

    # Plot 2: Size and operation of the FC
    ax2.plot(range(len(df_input)), [P_FC[t] / 1000 for t in range(len(df_input))], label='FC', color=palette[0])
    ax2.plot(range(len(df_input)), [S_FC / 1000] * len(df_input), label='FC size', linestyle=':', color=palette[2])
    # ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Power [kW]')
    # ax2.set_ylim([0, S_FC_max * 1.1 / 1000])
    ax2.set_ylim([0, S_FC / 1000 * 1.1])
    ax2.legend(loc='upper right')

    # Plot 3: Energy and TANK Size
    ax3.plot(range(len(df_input)), [E_TANK[t] / 3600000 for t in range(len(df_input))], label='TANK', color=palette[5])
    ax3.plot(range(len(df_input)), [S_TANK / 3600000] * len(df_input), label='TANK size', linestyle='-.', color=palette[6])
    ax3.set_xlabel('Time [h]')
    ax3.set_ylabel('Energy [kWh]')
    # ax3.set_ylim([0, (S_TANK_max * 1.1) / 3600000])
    ax3.set_ylim([0, S_TANK / 3600000 * 1.1])
    ax3.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

def plot_costs_and_prices(cost_inst, cost_op, cost, cost_maint, df_input):
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
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(15, 15))

    # Plot electricity prices
    ax1.plot(range(len(df_input)), df_input['price_Eur_MWh'], label='Import Price', color=palette[0])
    ax1.plot(range(len(df_input)), df_input['Price_DayAhed'], label='Export Price', color=palette[1])
    ax1.set_xlabel('Time [h]', fontsize=24)
    ax1.set_ylabel('Price [€/Wh]', fontsize=28)
    ax1.legend(fontsize=24)
    
    # Increase the font size for the tick labels
    ax1.tick_params(axis='both', which='major', labelsize=24)

    # Plot installation cost, operation cost, main cost, startup cost, and total cost
    bar_width = 0.4
    bar1 = ax2.bar(0, cost_inst / 1000, bar_width, label='Installation Cost', color=palette[2])
    bar2 = ax2.bar(1, cost_op   / 1000, bar_width, label='Operation Cost', color=palette[3])
    bar3 = ax2.bar(2, cost_maint / 1000, bar_width, label='Main Cost', color=palette[4])
    bar4 = ax2.bar(3, cost      / 1000, bar_width, label='Total Cost', color=palette[6])

    ax2.set_ylabel('Cost [k€/y]', fontsize=28)
    ax2.set_xticks([0, 1, 2, 3])
    ax2.set_xticklabels(['Installation Cost', 'Operation Cost', 'Maintenance Cost', 'Total Cost'], fontsize=24)
    ax2.tick_params(axis='both', which='major', labelsize=24)  # Set tick label font size for ax2

    plt.tight_layout()
    plt.show()

