# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:25:19 2024

@author: fism
"""
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

kWh2J = 3600*1000
J2kWh = 1 / (3600*1000)

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

#------------------------------------------------------------------------------
# Inputs and Demand plots
#------------------------------------------------------------------------------

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
    
# def plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours):
#     # Plot power generation, imported, and exported power
#     fig, ax = plt.subplots(1, 1, figsize=(20, 10))  # Use a single axis
#     ax.plot(range(len(df_input)), [P_PV[t] / 1000 for t in range(nHours)], label='PV Generation', color='orange')
#     ax.plot(range(len(df_input)), [P_imp[t] / 1000 for t in range(nHours)], label='Imported Power', color='blue')
#     ax.plot(range(len(df_input)), [-P_exp[t] / 1000 for t in range(nHours)], label='Exported Power', color='green')
    
#     ax.set_ylabel('Power [kW]')
#     ax.legend(loc='upper left')
#     ax.set_title('PV Generation, Imported, and Exported Power')
#     ax.set_xlabel('Time [h]')  # Set the x-label here

#     plt.show()
    
def plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours):
    df_plot = pd.DataFrame({
        'Time': range(nHours),
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
    
    return fig
    # fig.write_html('plot_power_generation.html', auto_open=True)

#------------------------------------------------------------------------------
# Component sizes plots
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------
# System operation plots
#------------------------------------------------------------------------------

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
    ax1.set_ylabel('Power [kW]', fontsize=24)
    # ax1.set_ylim([0, S_ELY_max * 1.1 / 1000])
    ax1.set_ylim([0, S_ELY / 1000 * 1.1])
    
    ax1.tick_params(axis='both', which='major', labelsize=24)
    ax1.legend(loc='upper right', fontsize=24)

    # Plot 2: Size and operation of the FC
    ax2.plot(range(len(df_input)), [P_FC[t] / 1000 for t in range(len(df_input))], label='FC', color=palette[0])
    ax2.plot(range(len(df_input)), [S_FC / 1000] * len(df_input), label='FC size', linestyle=':', color=palette[2])
    # ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Power [kW]', fontsize=24)
    # ax2.set_ylim([0, S_FC_max * 1.1 / 1000])
    ax2.set_ylim([0, S_FC / 1000 * 1.1])
    ax2.tick_params(axis='both', which='major', labelsize=24)
    ax2.legend(loc='upper right', fontsize=24)

    # Plot 3: Energy and TANK Size
    ax3.plot(range(len(df_input)), [E_TANK[t] / 3600000 for t in range(len(df_input))], label='TANK', color=palette[5])
    ax3.plot(range(len(df_input)), [S_TANK / 3600000] * len(df_input), label='TANK size', linestyle='-.', color=palette[6])
    ax3.set_xlabel('Time [h]', fontsize=24)
    ax3.set_ylabel('Energy [kWh]', fontsize=24)
    # ax3.set_ylim([0, (S_TANK_max * 1.1) / 3600000])
    # ax3.set_ylim([0, S_TANK / 3600000 * 1.1])
    ax3.set_ylim([0, 5000])
    ax3.tick_params(axis='both', which='major', labelsize=24)
    ax3.legend(loc='upper right', fontsize=24)

    plt.tight_layout()
    plt.show()

def plot_battery_operation(P_demand, P_imp, P_ch, P_ds, E_b, bat_params, C_b_kWh, nHours):
    
    # First set of tile plots
    SOC_min = bat_params['SOC_min']
    SOC_max = bat_params['SOC_max']
    power_supplied = [(P_ds[t] + P_imp[t] - P_ch[t])/1000 for t in range(nHours)]
    
    fig, axs = plt.subplots(3, 1, figsize=(12, 18))
    
    # Tile 1: Energy demand and sum of discharged battery power and import grid vs time
    axs[0].plot(P_demand/1000,       label='Energy Demand (kWh)', color='blue')
    axs[0].plot(power_supplied, label='Sum of Discharged and Imported Power minus charged (kWh)', color='red')
    #axs[0].set_title('Energy Demand and Supply Over Time - They should be equal; I guess I counted battery eff twice somewhere')
    axs[0].set_xlabel('Time (hours)', fontsize=24)
    axs[0].set_ylabel('Energy (kWh)', fontsize=24)
    axs[0].set_ylim([0, 800])
    
    axs[0].tick_params(axis='both', which='major', labelsize=24)
    axs[0].legend(fontsize=20)

    # Tile 2: Import power, charging power, and discharging power vs time
    axs[1].plot([P_imp[t] /1000 for t in range(nHours)], label='Import Power (kW)', color='green')
    axs[1].plot([P_ch[t]  /1000 for t in range(nHours)], label='Charging Power (kW)', color='orange')
    axs[1].plot([P_ds[t]  /1000 for t in range(nHours)], label='Discharging Power (kW)', color='purple')
    #axs[1].set_title('Power Flows Over Time')
    axs[1].set_xlabel('Time (hours)', fontsize=24)
    axs[1].set_ylabel('Power (kW)',   fontsize=24)
    axs[1].set_ylim([0, 400])
    
    axs[1].tick_params(axis='both', which='major', labelsize=24)
    axs[1].legend(fontsize=20)

    # Tile 3: Energy in the battery vs time
    axs[2].plot([E_b[t]/3600000 for t in range(nHours)], label='Battery Energy (kWh)', color='cyan')
    #axs[2].set_title('Battery State of Charge Over Time')
    axs[2].axhline(y=SOC_min * C_b_kWh, color='orange', linestyle='--', 
                   label=f'Min SOC ({SOC_min * 100}%) Capacity: {SOC_min * C_b_kWh:.2f} kWh')
    axs[2].axhline(y=SOC_max * C_b_kWh, color='green', linestyle='--', 
                   label=f'Max SOC ({SOC_max * 100}%) Capacity: {SOC_max * C_b_kWh:.2f} kWh')
    axs[2].axhline(y=C_b_kWh, color='red', linestyle='-', 
                   label=f'Selected Capacity: {C_b_kWh:.2f} kWh')
    axs[2].set_xlabel('Time (hours)', fontsize=24)
    axs[2].set_ylabel('Energy (kWh)', fontsize=24)
    axs[2].set_ylim([0, 500])
    
    axs[2].tick_params(axis='both', which='major', labelsize=24)
    axs[2].legend(fontsize=20)

    plt.tight_layout()
    plt.show()


def plot_WHR(results):
    # Convert the results dictionary into a pandas DataFrame
    df = pd.DataFrame(results)
    
    # Create a figure with 2 rows and 1 column
    fig = make_subplots(rows=2, cols=1,
                        subplot_titles=("Cooling Water Mass Flows", "Heat Recovery from the Electrolyzer and the Fuel Cell"),
                        shared_xaxes=True)  # Ensure the x-axis (time) is shared between subplots
    
    # Add Cooling Water Mass Flows to the first subplot
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_ELY'], name='m_cw_ELY', mode='lines'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_FC'], name='m_cw_FC', mode='lines'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_HT'], name='m_cw_HT', mode='lines'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_MT'], name='m_cw_MT', mode='lines'), row=1, col=1)

    # Add Heat Recovery to the second subplot
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_th_HT'], name='P_th_HT', stackgroup='one'), row=2, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_th_MT'], name='P_th_MT', stackgroup='one'), row=2, col=1)
    
    # Add 'P_ELY + P_FC' line plot to the second subplot
    df['P_ELY_plus_P_FC'] = df['P_ELY'] + df['P_FC']
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_ELY_plus_P_FC'], name='P_ELY + P_FC', mode='lines+markers'), row=2, col=1)

    # Update layout to include range slider and selector
    fig.update_layout(
        height=800, 
        title_text="WHR Analysis",
        xaxis2=dict(
            rangeslider=dict(
                visible=True
            ),
            type="date",
            rangeselector=dict(
                buttons=list([
                    dict(count=1, label="1m", step="month", stepmode="backward"),
                    dict(count=6, label="6m", step="month", stepmode="backward"),
                    dict(count=1, label="YTD", step="year", stepmode="todate"),
                    dict(count=1, label="1y", step="year", stepmode="backward"),
                    dict(step="all")
                ])
            )
        )
    )

    # Save the figure as HTML and automatically open it
    #fig.write_html("WHR_plot.html", auto_open=True)
    return fig

#------------------------------------------------------------------------------
# Cost distribution plots
#------------------------------------------------------------------------------

def plot_costs_and_prices(all_costs, df_input):
    """
    Parameters:
    - all_costs: dict, contains various cost components and revenues
    - df_input: DataFrame, input data including electricity prices
    """
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(15, 15))

    # Plot electricity prices
    ax1.plot(range(len(df_input)), df_input['price_Eur_MWh'], label='Import Price', color=palette[0])
    ax1.plot(range(len(df_input)), df_input['Price_DayAhed'], label='Export Price', color=palette[1])
    ax1.set_xlabel('Time [h]', fontsize=24)
    ax1.set_ylabel('Price [€/Wh]', fontsize=28)
    ax1.legend(fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)

    bar_width = 0.4

    # Define bars and labels for the legend
    bars_and_labels = [
        ("Exported Electricity", -all_costs["exported electricity"] / 1000, palette[2]),
        ("Revenues WHR", -all_costs["revenues WHR"] / 1000, palette[3]),
        ("Imported Electricity", all_costs["imported electricity"] / 1000, palette[4]),
        ("Grid Use Fees", all_costs["grid use fees"] / 1000, palette[5])
    ]

    # Plot operational cost bars with cumulative stacking
    cumulative_negative = cumulative_positive = 0
    for name, value, color in bars_and_labels:
        if value < 0:
            ax2.bar(1, value, bar_width, bottom=cumulative_negative, color=color, label=name)
            cumulative_negative += value
        else:
            ax2.bar(1, value, bar_width, bottom=cumulative_positive, color=color, label=name)
            cumulative_positive += value

    # Plot other costs without adding to legend (skip label argument)
    ax2.bar(0, all_costs["installation cost"] / 1000, bar_width, color='skyblue')
    ax2.bar(2, all_costs["maintenance cost"] / 1000, bar_width, color='lightgreen')
    ax2.bar(3, all_costs["TAC"] / 1000, bar_width, color='salmon')

    ax2.set_ylabel('Cost [k€/y]', fontsize=28)
    ax2.set_xticks([0, 1, 2, 3])
    ax2.set_xticklabels(['Installation Cost', 'Operational Cost', 'Maintenance Cost', 'Total Annual Cost'], fontsize=18)
    ax2.tick_params(axis='both', which='major', labelsize=24)

    # Customize legend to only include the four operational cost parameters
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[-4:], labels[-4:], fontsize=20, loc='upper left')  # Adjust to ensure only the last four items are included

    plt.tight_layout()
    plt.show()

# Function to calculate operational cost and create pie chart data
def costs_pie_chart(all_costs):
    operational_cost = all_costs['imported electricity'] + all_costs['grid use fees'] - all_costs['exported electricity'] - all_costs['revenues WHR']
    pie_chart_data = {
        'Categories': ['Installation Cost', 'Operational Cost', 'Maintenance Cost'],
        'Values': [all_costs['installation cost'], operational_cost, all_costs['maintenance cost']]
    }
    fig = px.pie(pie_chart_data, names='Categories', values='Values', title='Cost Overview')
    
    return fig
    # fig.write_html('costs_pie_chart.html', auto_open=True)
