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
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.signal import gaussian, convolve
from statsmodels.nonparametric.smoothers_lowess import lowess


kWh2J = 3600*1000
J2kWh = 1 / (3600*1000)

# Set Seaborn style
sns.set(style="whitegrid")
palette = sns.color_palette("husl", 8)  # You can choose a different palette if needed

# Increase font sizes for better readability on PowerPoint slides
plt.rcParams.update({'font.size': 24})  # Adjust global font size
plt.rcParams['axes.labelsize']  = 24
plt.rcParams['axes.titlesize']  = 24
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['legend.fontsize'] = 24

#------------------------------------------------------------------------------
# Inputs and Demand plots
#------------------------------------------------------------------------------
def exponential_moving_average(data, alpha=0.01):
    return pd.Series(data).ewm(alpha=alpha).mean().values

def gaussian_smoothing(data, window_size=50, std=10):
    gaussian_window = gaussian(window_size, std=std)
    gaussian_window /= gaussian_window.sum()
    smoothed = convolve(data, gaussian_window, mode='same')
    return smoothed

def lowess_smoothing(data, frac=0.003): # Adjust frac to change smoothness level
    return lowess(data, np.arange(len(data)), frac=frac, return_sorted=False)

def moving_average(data, window=150):
    return np.convolve(data, np.ones(window)/window, mode='valid')

def plot_heat_demand(df_heat_demand):
    # Set up the figure and axis
    fig, axes = plt.subplots(nrows=2, figsize=(20, 10), sharex=True)
    sns.set(style="whitegrid")
    palette = sns.color_palette("Set2")

    # Apply moving average smoothing
    smooth_zone1_35 = lowess_smoothing(df_heat_demand['Heating_Zone1_35degC_kW'].values)
    smooth_zone1_60 = lowess_smoothing(df_heat_demand['Heating_Zone1_60degC_kW'].values)
    smooth_zone2_35 = lowess_smoothing(df_heat_demand['Heating_Zone2_35degC_kW'].values)
    smooth_zone2_60 = lowess_smoothing(df_heat_demand['Heating_Zone2_60degC_kW'].values)
    smooth_zone3_35 = lowess_smoothing(df_heat_demand['Heating_Zone3_35degC_kW'].values)
    smooth_zone3_60 = lowess_smoothing(df_heat_demand['Heating_Zone3_60degC_kW'].values)
    
    # Adjust the day range for the smoothed data
    smoothed_hours = np.arange(len(smooth_zone1_35))

    # Stack plot for Space Heating Demand
    axes[0].stackplot(smoothed_hours, smooth_zone3_35, smooth_zone1_35, smooth_zone2_35,
                      labels=['Zone 3', 'Zone 1', 'Zone 2'],
                      colors=[palette[1], palette[0], palette[3]],
                      alpha=0.8, edgecolor='black')
    axes[0].set_title('(a)', fontsize=26)
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Demand (kWh)')
    axes[0].legend(loc='upper center', fontsize=24)
    axes[0].tick_params(axis='both', labelsize=22)
    axes[0].grid(False)

    # Stack plot for Domestic Hot Water Demand
    axes[1].stackplot(smoothed_hours, smooth_zone3_60, smooth_zone1_60, smooth_zone2_60,
                      labels=['Zone 3', 'Zone 1', 'Zone 2'],
                      colors=[palette[1], palette[0], palette[3]],
                      alpha=0.8, edgecolor='black')
    axes[1].set_title('(b)', fontsize=26)
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Demand (kWh)')
    axes[1].legend(loc='upper right', fontsize=24)
    axes[1].tick_params(axis='both', labelsize=22)
    axes[1].grid(False)

    # Automatically adjust the layout
    plt.tight_layout()

    # Display the plot
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

     

def pv_efficiency(df_pv):
    # Plot
    plt.figure(figsize=[10,8])
    pc = plt.scatter(df_pv['irradiance'], df_pv['eta_cell'], c=df_pv['T_amb'], cmap='jet')
    plt.colorbar(label='Temperature [°C]', ax=plt.gca()) # Adding a colorbar with a label
    pc.set_alpha(0.25) # Setting the alpha for the scatter plot
    plt.grid(alpha=0.5) # Adding grid to the plot with some transparency
    plt.ylim(bottom=0.16)  # # Setting the lower limit for y-axis, adjusted for clarity
    plt.xlabel('Irradiance [W/m²]')
    plt.ylabel('Relative efficiency [-]')
    plt.tight_layout()
    plt.show()
    

def plot_power_generation(results, df_input, nHours):

    # Prepare the data
    df = pd.DataFrame(results)
    df_plot = pd.DataFrame({
        'Time': range(nHours),
        'PV Generation [kW]': [p / 1000 for p in df['P_PV']],
        'Imported Power [kW]': [p / 1000 for p in df['P_imp']],
        'Exported Power [kW]': [-p / 1000 for p in df['P_exp']],
        'Demand [kW]': [p / 1000 for p in df['P_demand']]
    })

    # Define colors from Seaborn's Set2 palette

    colors = sns.color_palette("Set2", 5)  # Four colors for the four components

    # Calculate total sums in MWh for legend labels
    total_pv_gen_MWh   = df_plot['PV Generation [kW]'].sum() / 1000
    total_imported_MWh = df_plot['Imported Power [kW]'].sum() / 1000
    total_exported_MWh = df_plot['Exported Power [kW]'].sum() / 1000
    total_demand_MWh   = df_plot['Demand [kW]'].sum() / 1000

    labels = [
        f'PV Generation, total: {total_pv_gen_MWh:.2f} MWh',
        f'Imported Power, total: {total_imported_MWh:.2f} MWh',
        f'Exported Power, total: {total_exported_MWh:.2f} MWh',
        f'Demand, total: {total_demand_MWh:.2f} MWh'
    ]


    # Plot the data using Matplotlib and Seaborn colors
    plt.figure(figsize=(20, 10))  # Set the figure size

    # # Stack plot for PV Generation and Imported Power                        # Nice for weekly or monthly representations
    # plt.stackplot(df_plot['Time'],df_plot['Imported Power [kW]'], df_plot['PV Generation [kW]'], 
    #               labels=labels[:2],
    #               colors=colors[:2])
    
    # Fill the area for Exported Power
    plt.fill_between(df_plot['Time'], 0, df_plot['PV Generation [kW]'], label=labels[0], color=colors[1], alpha=1)
    plt.fill_between(df_plot['Time'], 0, df_plot['Imported Power [kW]'], label=labels[1], color=colors[0], alpha=1)
    plt.fill_between(df_plot['Time'], 0, df_plot['Exported Power [kW]'], label=labels[2], color=colors[2], alpha=1)

    # Plot Demand with a dotted line
    # plt.plot('Time', 'Demand [kW]', data=df_plot, color=colors[3], linestyle='--', linewidth=2)

    # Update labels for the remaining lines
    # plt.plot([], [], color=colors[1], label=labels[1])                  # For Imported Power
    # plt.plot([], [], color=colors[0], label=labels[0])                  # For PV Power
    # plt.plot([], [], color=colors[2], label=labels[2])                  # For Exported Power
    plt.plot([], [], color=colors[3], linestyle='--', label=labels[3])    # For Demand

    # Style the plot
    plt.title('PV Generation, Imported, Exported Power, and Demand Overview', fontsize=24)
    plt.xlabel('Time [h]', fontsize=24)
    plt.ylabel('Power [kW]', fontsize=24)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), shadow=True, ncol=2, fontsize=24)
    plt.grid(True)
    plt.tick_params(axis='both', which='major', labelsize=24)

    # Show or save the figure
    plt.show()

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
def plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input):  # without PWA
# def plot_HESS_results(P_PV, P_ELY, P_ELY_PWA, S_ELY, S_ELY_max, P_FC, P_FC_in, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input): #if PWA
    """
    Plot main results including the size and operation of the electrolyzer (ELY),
    fuel cell (FC), and energy and TANK size over time.

    Parameters:
    - P_PV: list, hourly PV power generation values in [W]
    - P_ELY: list, hourly electrolyzer power values in [W]
    - S_ELY: float, electrolyzer size in [W]
    - P_FC: list, hourly fuel cell power values in [W]
    - S_FC: float, fuel cell size in [W]
    - E_TANK: list, hourly tank energy values in [J]
    - S_TANK: float, tank size
    - S_TANK_max: float, maximal tank size
    - df_input: DataFrame, input data
    """
    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(20, 10))

    # Plot 1: Size and operation of the ELY
    ax1.plot(range(len(df_input)), [P_ELY[t] / 1000 for t in range(len(df_input))], label='P_ELY_in', color=palette[1])
    # ax1.plot(range(len(df_input)), [P_ELY_PWA[t] / 1000 for t in range(len(df_input))], label='P_ELY_out', color=palette[4])
    ax1.plot(range(len(df_input)), [S_ELY / 1000] * len(df_input), 
             label=f'ELY size: {S_ELY / 1000:.2f} kW', linestyle='--', color=palette[3])
    # ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Power [kW]', fontsize=24)
    # ax1.set_ylim([0, S_ELY_max * 1.1 / 1000])
    ax1.set_ylim([0, S_ELY / 1000 * 1.1])
    
    ax1.tick_params(axis='both', which='major', labelsize=24)
    ax1.legend(loc='upper right', fontsize=24)

    # Plot 2: Size and operation of the FC
    # ax2.plot(range(len(df_input)), [P_FC_in[t] / 1000 for t in range(len(df_input))], label='P_FC_in', color=palette[0])
    ax2.plot(range(len(df_input)), [P_FC[t]    / 1000 for t in range(len(df_input))], label='P_FC_out', color=palette[6])
    ax2.plot(range(len(df_input)), [S_FC / 1000] * len(df_input), 
             label=f'FC size: {S_FC/ 1000:.2f} kW', linestyle=':', color=palette[2])
    # ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Power [kW]', fontsize=24)
    # ax2.set_ylim([0, S_FC_max * 1.1 / 1000])
    ax2.set_ylim([0, S_FC / 1000 * 1.1])
    ax2.tick_params(axis='both', which='major', labelsize=24)
    ax2.legend(loc='upper right', fontsize=24)

    # Plot 3: Energy and TANK Size
    ax3.fill_between(range(len(df_input)), 0, [E_TANK[t] / 3600000 for t in range(len(df_input))], label='TANK', color=palette[5])
    # ax3.plot(range(len(df_input)), [E_TANK[t] / 3600000 for t in range(len(df_input))], label='TANK', color=palette[5])
    ax3.plot(range(len(df_input)), [S_TANK / 3600000] * len(df_input), 
             label=f'H2-TANK size: {S_TANK/ (39.39*3600000):.2f} kgH2', linestyle='-.', color=palette[6])  
    ax3.set_xlabel('Time [h]', fontsize=24)
    ax3.set_ylabel('Energy [kWh]', fontsize=24)
    # ax3.set_ylim([0, (S_TANK_max * 1.1) / 3600000])
    # ax3.set_ylim([0, S_TANK / 3600000 * 1.1])
    ax3.set_ylim([0, 35000])
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
    
    # Divide specified columns by 1000 to convert values from Watts to Kilowatts
    df['P_th_HT'] = df['P_th_HT'] / 1000
    df['P_th_MT'] = df['P_th_MT'] / 1000
    df['P_ELY']   = df['P_ELY'] / 1000
    df['P_FC']   = df['P_FC'] / 1000
    
    # Create a figure with 2 rows and 1 column, adjusting subplot positions
    fig = make_subplots(rows=2, cols=1,
                        subplot_titles=("Heat Recovery", "Cooling Water Mass Flows"),
                        shared_xaxes=True)  # Ensure the x-axis (time) is shared between subplots
    
    # First subplot for Heat Recovery in Kilowatt-hours
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_th_HT'], name='P_th_HT', stackgroup='one'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_th_MT'], name='P_th_MT', stackgroup='one'), row=1, col=1)
    df['P_ELY_plus_P_FC'] = df['P_ELY'] + df['P_FC']  # Sum already in Kilowatts after division
    fig.add_trace(go.Scatter(x=df['ts'], y=df['P_ELY_plus_P_FC'], name='P_ELY + P_FC', mode='lines+markers'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['z3_35degC_kWh'], name='Heat Demand Zone 3 @ 35°C', mode='lines', line=dict(dash='dot')), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['z3_60degC_kWh'], name='Heat Demand Zone 3 @ 60°C', mode='lines', line=dict(dash='dot')), row=1, col=1)

    # Second subplot for Cooling Water Mass Flows in Kilograms per second
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_ELY'], name='m_cw_ELY', mode='lines'), row=2, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_FC'], name='m_cw_FC', mode='lines'), row=2, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_HT'], name='m_cw_HT', mode='lines'), row=2, col=1)
    fig.add_trace(go.Scatter(x=df['ts'], y=df['m_cw_MT'], name='m_cw_MT', mode='lines'), row=2, col=1)

    # Update layout to include range slider and selector, specify units on the Y-axes
    fig.update_layout(
        height=800, 
        title_text="WHR Analysis",
        xaxis=dict(
            rangeslider=dict(
                visible=True
            ),
            type="date",
            rangeselector=dict(
                buttons=[
                    dict(count=1, label="1y", step="year", stepmode="todate"),
                    dict(count=1, label="1m", step="month", stepmode="todate"),
                    dict(count=1, label="1d", step="day", stepmode="todate"),
                ]
            )
        ),
        xaxis2=dict(
            type="date"
        )
    )
    
    # Update y-axes titles to indicate units
    fig.update_yaxes(title_text="Kilowatt-hours (kWh)", row=1, col=1)
    fig.update_yaxes(title_text="Kilograms per second (kg/s)", row=2, col=1)

    # Save the figure as HTML and automatically open it
    #fig.write_html("WHR_plot.html", auto_open=True)
    return fig

#------------------------------------------------------------------------------
# Cost distribution plots
#------------------------------------------------------------------------------

def plot_costs_and_prices(all_costs, df_input, electricity_price_imp, electricity_price_exp):
    """
    Parameters:
    - all_costs: dict, contains various cost components and revenues
    - df_input: DataFrame, input data including electricity prices
    """
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(15, 15))

    # Plot electricity prices
    ax1.plot(range(len(df_input)), df_input['price_Eur_MWh'], label='Spot Market Import Price in €', color=palette[0])
    ax1.plot(range(len(df_input)), electricity_price_imp, label='BKW import price in CHF', color=palette[2])
    ax1.plot(range(len(df_input)), df_input['Price_DayAhed'], label='Day Ahead Export Price in €', color=palette[1])
    ax1.plot(range(len(df_input)), electricity_price_exp, label='BKW export price in CHF', color=palette[4])
    ax1.set_xlabel('Time [h]', fontsize=24)
    ax1.set_ylabel('Price per MWh', fontsize=28)
    ax1.legend(fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)

    bar_width = 0.4

    # Define bars and labels for the legend
    bars_and_labels = [
        ("Exported Electricity", -all_costs["operational cost"]["exported electricity"] / 1000, palette[2]),
        ("Revenues WHR", -all_costs["revenues WHR"] / 1000, palette[6]),
        ("Imported Electricity", all_costs["operational cost"]["imported electricity"] / 1000, palette[4]),
        ("Grid Use Fees", all_costs["operational cost"]["grid use fees"]["grid_kWh"] / 1000, palette[5]),
        ("Taxes & Levies", all_costs["operational cost"]["taxes & levies"] / 1000, palette[7])
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
    ax2.bar(0, all_costs["installation cost"]["total"] / 1000, bar_width, color='skyblue')
    ax2.bar(2, all_costs["maintenance cost"]["total"] / 1000, bar_width, color='lightgreen')
    ax2.bar(3, all_costs["TAC"] / 1000, bar_width, color='salmon')

    ax2.set_ylabel('Cost [k€/y]', fontsize=28)
    ax2.set_xticks([0, 1, 2, 3])
    ax2.set_xticklabels(['Installation Cost', 'Operational Cost', 'Maintenance Cost', 'Total Annual Cost'], fontsize=18)
    ax2.tick_params(axis='both', which='major', labelsize=24)

    # Customize legend to only include the four operational cost parameters
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[-5:], labels[-5:], fontsize=20, loc='upper left')  # Adjust to ensure only the last four items are included

    plt.tight_layout()
    plt.show()

# Function to calculate operational cost and create pie chart data
def costs_pie_chart(all_costs):
    operational_cost = all_costs['imported electricity'] + all_costs['grid use fees'] + all_costs['taxes & levies']  - all_costs['exported electricity'] - all_costs['revenues WHR']
    pie_chart_data = {
        'Categories': ['Installation Cost', 'Operational Cost', 'Maintenance Cost'],
        'Values': [all_costs['installation cost'], operational_cost, all_costs['maintenance cost']]
    }
    fig = px.pie(pie_chart_data, names='Categories', values='Values', title='Cost Overview')
    
    return fig
    # fig.write_html('costs_pie_chart.html', auto_open=True)
    
#------------------------------------------------------------------------------
# ELY PWA Plots
#------------------------------------------------------------------------------
def plot_efficiencies(P_ELY, P_ELY_PWA, S_ELY, nHours, x_bp_val, y_bp_val_scaled, Vdot_FC_H2, i_FC, P_FC_in, P_FC, S_FC):
    sns.set_theme(style="whitegrid")
    
    # Calculating & sorting input & output powers and efficiencies for the Electrolyser
    inputPower  = sorted([P_ELY[t] / S_ELY for t in range(nHours)])
    outputPower = sorted([P_ELY_PWA[t] / S_ELY for t in range(nHours)])
    eta_ELY     = [(outp / inp) * 100 if inp != 0 else 0 for inp, outp in zip(inputPower, outputPower)]
    eta_bp      = [(y / x) * 100 for x, y in zip(x_bp_val, y_bp_val_scaled)] 

    # Calculate efficiency in % for the fuel cell
    P_FC_in_kW  = sorted([P_FC_in[t] / (S_FC) for t in range(nHours)])
    P_FC_out_kW = sorted([P_FC[t]    / (S_FC) for t in range(nHours)])
    eta_FC      = [(out / inp) * 100 if inp != 0 else 0 for inp, out in zip(P_FC_in_kW, P_FC_out_kW)]

    # Create a figure with a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 17))

    # First subplot: Electrolyzer Input/Output Power Curve
    norm = mcolors.Normalize(vmin=57, vmax=67)
    sc1 = axs[0, 0].scatter(inputPower, outputPower, c=eta_ELY, cmap='viridis', norm=norm, marker='o')
    cbar1 = fig.colorbar(sc1, ax=axs[0, 0])
    cbar1.set_label('Efficiency [%]', fontsize=20)
    axs[0, 0].set_xlabel('Input Power', fontsize=20)
    axs[0, 0].set_ylabel('Output Power', fontsize=20)
    axs[0, 0].set_title('ELY Normalized Input/Output Power Curve', fontsize=20)
    axs[0, 0].grid(True)
    axs[0, 0].tick_params(axis='both', labelsize=20)

    # Second subplot: Electrolyzer Efficiency vs. Input Power
    axs[0, 1].scatter(inputPower, eta_ELY, color='blue', marker='o')
    axs[0, 1].plot(x_bp_val, eta_bp, marker='o', color='red', label='eta_bp')  # Red color for distinction
    axs[0, 1].set_xlabel('Input Power', fontsize=20)
    axs[0, 1].set_ylabel('Efficiency [%]', fontsize=20)
    axs[0, 1].set_xlim(0, 1.2)
    axs[0, 1].set_ylim(40, 70)
    axs[0, 1].set_title('Electrolyser Efficiency', fontsize=20)
    axs[0, 1].grid(True)
    axs[0, 1].tick_params(axis='both', labelsize=20)
    axs[0, 1].legend()

    # Third subplot: P_FC_in_kW vs P_FC_out_kW for the fuel cell (with colorbar)
    norm = mcolors.Normalize(vmin=45, vmax=55)
    sc2   = axs[1, 0].scatter(P_FC_in_kW, P_FC_out_kW, c=eta_FC, cmap='plasma', norm=norm)
    cbar2 = fig.colorbar(sc2, ax=axs[1, 0])
    cbar2.set_label('Efficiency [%]', fontsize=20)
    axs[1, 0].set_xlabel('Input Power (kW)', fontsize=20)
    axs[1, 0].set_ylabel('Output Power (kW)', fontsize=20)
    axs[1, 0].tick_params(axis='both', labelsize=20)
    axs[1, 0].set_title('FC Input/Output Power Curve', fontsize=20)

    # Fourth subplot: i_FC vs Vdot_FC_H2 for the fuel cell (no colorbar)
    axs[1, 1].scatter(Vdot_FC_H2, i_FC, c='blue')  # Using a solid color for simplicity
    axs[1, 1].set_xlabel('Vdot_FC_H2', fontsize=20)
    axs[1, 1].set_ylabel('i_FC', fontsize=20)
    axs[1, 1].tick_params(axis='both', labelsize=20)
    axs[1, 1].set_title('Fuel cell current to H2 volume flow', fontsize=20)
