# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:54:21 2024

@author: fism
"""

import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import numpy as np

# # Path to the uploaded Excel file
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-07_09-58-27 - current - SS75.xlsx'

# # Load the data from the specified sheet
# df_optimization_results = pd.read_excel(file_path, sheet_name='Optimization Results')

# # Assuming the dataset is hourly and covers a full year (8760 hours)
# # We need to create a datetime index for 2019
# date_range = pd.date_range(start='2019-01-01', periods=8760, freq='H')
# df_optimization_results['timestamp'] = date_range
# df_optimization_results.set_index('timestamp', inplace=True)

# # Convert power values from Watts to kilowatts
# df_optimization_results[['P_demand', 'P_imp', 'P_PV', 'P_exp', 'P_FC']] /= 1000

# # Selecting specific weeks
# weeks = {
#     "February Week": ('2019-02-06', '2019-02-08'),
#     "May Week": ('2019-05-08', '2019-05-10'),
#     "August Week": ('2019-08-07', '2019-08-09'),
#     "November Week": ('2019-11-06', '2019-11-08')
# }

# # Set up the matplotlib figure
# plt.figure(figsize=(22, 14))
# sns.set(style="white")  # Changed from 'whitegrid' to 'white' to remove grid lines

# # Loop through each specified week and plot
# for i, (title, (start, end)) in enumerate(weeks.items(), 1):
#     plt.subplot(2, 2, i)
#     week_data = df_optimization_results.loc[start:end]
#     sns.lineplot(data=week_data, x=week_data.index, y='P_demand', label='Demand', color='blue')
#     plt.fill_between(week_data.index, 0, week_data['P_imp'], color='red', alpha=0.4, label='Import')
#     plt.fill_between(week_data.index, week_data['P_imp'], week_data['P_imp'] + week_data['P_PV'], color='green', alpha=0.4, label='PV')
#     plt.fill_between(week_data.index, week_data['P_imp'] + week_data['P_PV'], week_data['P_imp'] + week_data['P_PV'] + week_data['P_FC'], color='orange', alpha=0.7, label='Fuel Cell')
#     plt.fill_between(week_data.index, 0, -week_data['P_exp'], color='purple', alpha=0.4, label='Export')

#     plt.xlabel('Time', fontsize=24)
#     plt.ylabel('Power (kW)', fontsize=24)
#     plt.title(title, fontsize=24)
#     plt.xticks(rotation=45, fontsize=20)
#     plt.yticks(fontsize=20)

# # Place the legend centrally below the plots
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), fancybox=True, shadow=True, ncol=5, fontsize=18)

# plt.tight_layout()
# plt.savefig('SS75-current.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# PLOT HEAT DEMAND
#------------------------------------------------------------------------------
# Data
data = {
    "Scenario": ["Grid Only", "Base Case", "SS75%", "100kW BAT", "1MW BAT"],
    "Maximal Import Peak in kW": [613.57, 562.67, 410.73, 391.53, 214.75]
}

# Create DataFrame
df = pd.DataFrame(data)

# Calculate percentage decrease
base_value = df['Maximal Import Peak in kW'][0]  # First bar value
df['Percentage Decrease'] = - ((base_value - df['Maximal Import Peak in kW']) / base_value) * 100


# Plotting
plt.figure(figsize=(7, 6))
barplot = sns.barplot(x="Scenario", y="Maximal Import Peak in kW", data=df, color='skyblue')

# Trend line
# Calculate means for each scenario (as they are unique here)
x = np.arange(len(df['Scenario']))  # the label locations
y = df['Maximal Import Peak in kW']
# z = np.polyfit(x, y, 1)  # Linear fit
# p = np.poly1d(z)
# plt.plot(x, p(x), "r--")  # Plotting the trend line

# Adding percentage decrease annotations
for i in range(1, len(df)):
    plt.text(x=i, y=df['Maximal Import Peak in kW'][i]+10, # Position the text slightly above the bar
              s= f"{df['Percentage Decrease'][i]:.2f}%", 
              color='red', 
              ha='center', fontsize=18)

# Adding labels
plt.xlabel('Scenario', fontsize=20)
plt.ylabel('Peak in kW', fontsize=20)
plt.title('Maximal Import Peak by Scenario', fontsize=20)
plt.xticks(x, df['Scenario'], fontsize=16)  # Ensure scenario names are shown
plt.yticks(fontsize=20)

plt.tight_layout()
plt.savefig('peak-shaving.png', bbox_inches='tight', dpi=150)
plt.show()

#------------------------------------------------------------------------------
# PLOT HEAT DEMAND
#------------------------------------------------------------------------------
# # Path to the uploaded Excel file
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-06_22-10-13 - current - PV5000 - Base Case Scenario.xlsx'

# # Load the data from the Excel file
# df_optimization_results = pd.read_excel(file_path, sheet_name='Optimization Results')

# # Assuming the dataset is hourly and covers a full year (8760 hours)
# # We need to create a datetime index for the year 2019 (or the relevant year)
# date_range = pd.date_range(start='2019-01-01', periods=8760, freq='H')
# df_optimization_results['timestamp'] = date_range
# df_optimization_results.set_index('timestamp', inplace=True)

# df_optimization_results[['z1_35degC_kWh', 'z2_35degC_kWh', 'z3_35degC_kWh', 'z1_60degC_kWh', 'z2_60degC_kWh', 'z3_60degC_kWh']] /= 1000

# # Splitting data by temperature categories
# df_35degC = df_optimization_results[['z1_35degC_kWh', 'z2_35degC_kWh', 'z3_35degC_kWh']] 
# df_60degC = df_optimization_results[['z1_60degC_kWh', 'z2_60degC_kWh', 'z3_60degC_kWh']] 

# # Resample data to daily sums
# daily_data_35degC = df_optimization_results[['z1_35degC_kWh', 'z2_35degC_kWh', 'z3_35degC_kWh']].resample('D').sum()
# daily_data_60degC = df_optimization_results[['z1_60degC_kWh', 'z2_60degC_kWh', 'z3_60degC_kWh']].resample('D').sum()

# # Assuming df_optimization_results is loaded as shown in your example

# # Set seaborn style
# sns.set(style="white")
# sns.set_context("talk")  # Scale fonts up

# # Define specific colors from the 'husl' palette
# husl_palette = sns.color_palette("husl", 8)  # Generate enough colors
# custom_colors = [husl_palette[4], husl_palette[1], husl_palette[0]]  # Colors at indices 5, 0, and 4

# # Creating the subplots
# fig, axs = plt.subplots(2, 1, figsize=(16, 14), sharex=True)

# # Function to plot stacked data using fill_between
# def plot_stacked(df, ax, title, show_legend=False):
#     # Explicitly naming zones according to the data columns
#     zone_labels = ['Zone 1', 'Zone 2', 'Zone 3']
#     colors = custom_colors
#     bottom = None
#     for column, label, color in zip(df.columns, zone_labels, colors):
#         data = df[column]
#         ax.fill_between(df.index, bottom if bottom is not None else 0, data + (bottom if bottom is not None else 0), alpha=0.9, label=label, color=color)
#         bottom = data if bottom is None else bottom + data

#     ax.set_title(title, fontsize=30)
#     ax.set_ylabel('Demand (kWh)', fontsize=28)
#     ax.set_xlabel('Time', fontsize=28)
#     if show_legend:
#         ax.legend(title='Zone', loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, title_fontsize='28', fontsize=28)
#         # Explicitly set the fontsize for y-axis ticks
#     ax.tick_params(axis='y', labelsize=28)  # Ensure y-axis ticks are the correct size
# # Plotting
# plot_stacked(daily_data_35degC, axs[0], 'Space Heating')
# plot_stacked(daily_data_60degC, axs[1], 'Domestic Hot Water', show_legend=True)

# plt.xticks(fontsize=28)
# plt.yticks(fontsize=28)
# plt.tight_layout() 
# plt.savefig('heat-demand.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# BASE CASE SCENARIO
#------------------------------------------------------------------------------
"""
Using data from Grid only, PV2500, Base Case ans PV7500
Plotting TAC, Installation, Operation & Maintenance costs over Area PV

"""
# # Data from the user's table
# pv_areas           = np.array([0, 2500, 5000, 7500])             # Area_PV_max in m2
# tac_values         = [285.57, 205.17, 166.35, 140.70]  # TAC (CHF/y)
# installation_costs = [0.00, 28.29, 56.58, 84.87]
# operational_costs  = [285.57, 175.47, 106.95, 51.59]   # Includes operational cost_total
# export_costs       = [0.00, -13.64, -48.73, -89.78]    # Negative values as they are credits (exports)
# maintenance_costs  = [0.00, 1.41, 2.83, 4.24]

# # Set the color palette
# colors = sns.color_palette('Set2')

# plt.figure(figsize=(10, 6))

# # Plotting each line with a thicker line width and different markers
# plt.plot(pv_areas, tac_values,         marker='s', linestyle='-', color=colors[0], label='TAC', linewidth=3, markersize=12)
# plt.plot(pv_areas, installation_costs, marker='^', linestyle='-', color=colors[1], label='Installation', linewidth=3, markersize=12)
# plt.plot(pv_areas, operational_costs,  marker='o', linestyle='-', color=colors[2], label='Operations (incl. Exports)', linewidth=3, markersize=12)
# plt.plot(pv_areas, export_costs,       marker='d', linestyle=':', color=colors[3], label='Exports', linewidth=3, markersize=12)  # Dotted line
# plt.plot(pv_areas, maintenance_costs,  marker='p', linestyle='-', color=colors[4], label='Maintenance', linewidth=3, markersize=12)

# # Adjusting font sizes for the labels, title, and ticks
# plt.xlabel('PV Area in m²', fontsize=20)
# plt.ylabel('Cost in kCHF/Year', fontsize=20)
# # plt.title('Cost Breakdown vs. PV Area', fontsize=14)
# plt.xticks(np.arange(0, 8000, step=2500), fontsize=18)
# plt.yticks(fontsize=18)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=18)

# # Adjusting the legend with a larger font size
# plt.grid(False)

# # Show the plot
# plt.tight_layout() 
# plt.savefig('base-case.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# Path to the uploaded Excel file
file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-07_09-58-27 - current - SS75.xlsx'
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-06_22-10-13 - current - PV5000 - Base Case Scenario.xlsx'
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-07_03-18-17 - current - with Battery.xlsx'
# Load the data from the Excel file
df_optimization_results = pd.read_excel(file_path, sheet_name='Optimization Results')

# Assuming the dataset is hourly and covers a full year (8760 hours)
# We need to create a datetime index for the year 2019 (or the relevant year)
date_range = pd.date_range(start='2019-01-01', periods=8760, freq='H')
df_optimization_results['timestamp'] = date_range
df_optimization_results.set_index('timestamp', inplace=True)

# Convert power values from Watts to kilowatts
# df_optimization_results[['P_demand', 'P_imp', 'P_PV', 'P_exp', 'P_FC', 'P_ch', 'P_ds']] /= 1000
df_optimization_results[['P_demand', 'P_imp', 'P_PV', 'P_exp', 'P_FC']] /= 1000

# Resample data to daily sums
# daily_data = df_optimization_results[['P_demand', 'P_imp', 'P_PV', 'P_exp', 'P_FC', 'P_ch', 'P_ds']].resample('D').sum()
daily_data = df_optimization_results[['P_demand', 'P_imp', 'P_PV', 'P_exp', 'P_FC']].resample('D').sum()

# #------------------------------------------------------------------------------
# #FOR 1 YEAR 
# #Set up the matplotlib figure
# plt.figure(figsize=(20, 8))
# sns.set(style="white")

# # Plotting the whole year data
# sns.lineplot(data=daily_data, x=daily_data.index, y='P_demand', label='Demand', color='blue')
# plt.fill_between(daily_data.index, 0, daily_data['P_imp'], color='red', alpha=0.4, label='Import')
# plt.fill_between(daily_data.index, daily_data['P_imp'], daily_data['P_imp'] + daily_data['P_PV'], color='green', alpha=0.4, label='PV')
# plt.fill_between(daily_data.index, 0, -daily_data['P_exp'], color='purple', alpha=0.4, label='Export')

# plt.xlabel('Date', fontsize=24)
# plt.ylabel('Power (kW)', fontsize=24)
# plt.title('Yearly Power Demand and Supply for the base case scenario 2019', fontsize=26)
# plt.xticks(fontsize=24)
# plt.yticks(fontsize=24)

# # Legend configuration
# plt.legend(loc='lower right', fontsize=22)

# plt.tight_layout()
# plt.savefig('Yearly-power-consumption-2019.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
#FOR MONTHLY TIMELINE
# All months

# # Determine uniform y-axis limits
# y_max = daily_data.max().max()
# y_min = -(daily_data['P_exp'].abs().max() * 1.1)  # Calculate y_min based on maximum of P_exp scaled up by 10%

# # Set up the matplotlib figure
# plt.figure(figsize=(24, 28))
# sns.set(style="white")

# # Month names for titles
# months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

# # Loop through each month and plot
# for month in range(1, 13):
#     plt.subplot(4, 3, month)
#     month_data = daily_data[daily_data.index.month == month]
#     sns.lineplot(data=month_data, x=month_data.index, y='P_demand', label='Demand', color='blue')
#     plt.fill_between(month_data.index, 0, month_data['P_imp'], color='red', alpha=0.4, label='Import')
#     plt.fill_between(month_data.index, month_data['P_imp'], month_data['P_imp'] + month_data['P_PV'], color='green', alpha=0.4, label='PV')
#     plt.fill_between(month_data.index, 0, -month_data['P_exp'], color='purple', alpha=0.4, label='Export')

#     # plt.xlabel('Date', fontsize=28)
#     plt.ylabel('Power (kW)', fontsize=28)
#     plt.title(months[month-1], fontsize=28)
#     plt.xticks(rotation=45, fontsize=22)
#     plt.yticks(fontsize=24)
#     plt.ylim(y_min, y_max)  # Apply the uniform y-axis limits

# # Place the legend centrally below the plots
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20), fancybox=True, shadow=True, ncol=4, fontsize=14)

# plt.tight_layout()
# plt.savefig('Yearly-power-consumption-by-month.png', bbox_inches='tight', dpi=150)
# plt.show()
#------------------------------------------------------------------------------
# #January and July
# # Determine uniform y-axis limits
# y_max = daily_data.max().max()
# y_min = -(daily_data['P_exp'].abs().max() * 1.1)  # Calculate y_min based on maximum of P_exp scaled up by 10%

# # Set up the matplotlib figure, two plots in one row
# plt.figure(figsize=(12, 6))  # Adjusted for a more suitable size for a single row
# sns.set(style="white")

# # Months to plot
# months_to_plot = [1, 7]  # January is 1, July is 7
# month_names = ['January', 'July']

# # Plotting January and July
# for i, month in enumerate(months_to_plot):
#     plt.subplot(1, 2, i + 1)  # One row, two columns, plot position i+1
#     month_data = daily_data[daily_data.index.month == month]
#     sns.lineplot(data=month_data, x=month_data.index, y='P_demand', color='blue')
#     plt.fill_between(month_data.index, 0, month_data['P_imp'], color='red', alpha=0.4, label='Import')
#     plt.fill_between(month_data.index, month_data['P_imp'], month_data['P_imp'] + month_data['P_PV'], color='green', alpha=0.4, label='PV')
#     plt.fill_between(month_data.index, 0, -month_data['P_exp'], color='purple', alpha=0.4, label='Export')

#     plt.ylabel('Power (kW)', fontsize=20)
#     plt.xlabel('Date', fontsize=20)
#     plt.title(month_names[i], fontsize=20)
#     plt.xticks(rotation=45, fontsize=16)
#     plt.yticks(fontsize=20)
#     plt.ylim(y_min, y_max)  # Apply the uniform y-axis limits

# # Place the legend at the bottom of the figure
# # plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=4, fontsize=20)

# plt.tight_layout()
# plt.savefig('Jan-Jul-power-consumption.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
#FOR 4 MONTH with H2
# # Selecting specific weeks
# weeks = {
#     "February Week": ('2019-02-06', '2019-02-08'),
#     # "May Week": ('2019-05-08', '2019-05-10'),
#     "August Week": ('2019-08-07', '2019-08-09'),
#     # "November Week": ('2019-11-06', '2019-11-08')
# }

# # Determine uniform y-axis limits
# y_max = 650
# y_min = -330

# # Set up the matplotlib figure
# plt.figure(figsize=(18, 14))
# sns.set(style="white")  # Changed from 'whitegrid' to 'white' to remove grid lines

# # Loop through each specified week and plot
# for i, (title, (start, end)) in enumerate(weeks.items(), 1):
#     plt.subplot(2, 2, i)
#     week_data = df_optimization_results.loc[start:end]
#     sns.lineplot(data=week_data, x=week_data.index, y='P_demand', label='Demand', color='blue')
#     plt.fill_between(week_data.index, 0, week_data['P_imp'], color='red', alpha=0.4, label='Import')
#     plt.fill_between(week_data.index, week_data['P_imp'], week_data['P_imp'] + week_data['P_PV'], color='green', alpha=0.4, label='PV')
#     plt.fill_between(week_data.index, week_data['P_imp'] + week_data['P_PV'], week_data['P_imp'] + week_data['P_PV'] + week_data['P_FC'], color='orange', alpha=0.7, label='Fuel Cell')
#     # plt.fill_between(week_data.index, week_data['P_imp'] + week_data['P_PV'], week_data['P_imp'] + week_data['P_PV'] + week_data['P_ds'], color='orange', alpha=0.7, label='Battery')
#     plt.fill_between(week_data.index, 0, -week_data['P_exp'], color='purple', alpha=0.4, label='Export')

#     plt.xlabel('Time', fontsize=24)
#     plt.ylabel('Power (kW)', fontsize=24)
#     plt.title(title, fontsize=24)
#     plt.xticks(rotation=45, fontsize=20)
#     plt.yticks(fontsize=20)
#     plt.ylim(y_min, y_max)  # Apply the uniform y-axis limits

# # Place the legend centrally below the plots
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.30), fancybox=True, shadow=True, ncol=5, fontsize=18)

# plt.tight_layout()
# plt.savefig('SS75-current.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# Self-Sufficiency
#------------------------------------------------------------------------------
# Data
data = {
    "Self-Sufficiency in %": [0, 25, 50, 75],
    "low efficiencies high costs": [118.12, 118.12, 118.13, 222.51],
    "improved efficiencies low costs": [118.12, 118.12, 118.12, 145.50],
    # "low efficiencies high costs w/o Fees for import peaks": [63.07, 63.07, 63.07, 182.03],
    "low efficiencies high costs w/o Grid Charges": [53.29, 53.29, 53.29, 176.45],
}

# Create DataFrame
df = pd.DataFrame(data)

# Additional points
additional_points = pd.DataFrame({
    "Self-Sufficiency in %": [75, 75],
    "Value": [107.51, 224.44],
    "Scenario": ["1MW Battery", "100kW Battery"]
})

# Seaborn styling
sns.set(style="white", palette="Set2")

# Create plot
plt.figure(figsize=(8, 8))
markers = ['v', '^', 's', 'P']#, '*']  # Different markers for each line

for i, column in enumerate(df.columns[1:]):
    sns.lineplot(x=df["Self-Sufficiency in %"], y=df[column], marker=markers[i], markersize=12, linewidth=4, label=column) # ,label=column

# Custom colors for each battery capacity point
colors = ['red', 'orange']
labels = ['1MWh Battery', 'HESS + 100kWh Battery']

# Plot additional points with individual labels and colors
for idx, row in additional_points.iterrows():
    plt.scatter(row["Self-Sufficiency in %"], row["Value"], color=colors[idx], s=200, label=labels[idx])

# plt.title('LCOE vs Self-Sufficiency', fontsize=28)
plt.xlabel('Self-Sufficiency in %', fontsize=24)
plt.ylabel('LCOE in CHF/MWh', fontsize=24)
plt.xticks([0, 25, 50, 75], fontsize=24)  # Set specific steps on x-axis
plt.yticks(fontsize=24)
plt.ylim(0, 250)  # Set y-axis minimum limit to 0
plt.legend(title='Scenario', title_fontsize='24', fontsize='24', loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig('Self-Sufficiency.png', bbox_inches='tight', dpi=150)
plt.show()


#------------------------------------------------------------------------------
# IMPACT OF GRID CHARGES
#------------------------------------------------------------------------------

# # Creating a DataFrame with the provided data
# data_updated = {
#     'Installation': [203, 198],
#     'Imported Electricity': [24, 25],
#     'Exported Electricity': [-10, -10],  # Negative values
#     'Grid Tariff': [8, 0],
#     'Fees for Peak Imports': [46, 0],
#     'Taxes': [13, 13],
#     'Maintenance': [12, 11],
#     'Revenues WHR': [-16, -15]  # Negative values
# }

# # Convert the dictionary to a DataFrame
# df_updated = pd.DataFrame(data_updated, index=['SS75', 'w/o Grid Charges'])

# # Plotting the updated stacked bar chart
# fig, ax = plt.subplots(figsize=(10, 10))
# # Setting the color palette to 'husl'
# sns.set(style="white", palette="Paired")

# df_updated.plot(kind='bar', stacked=True, ax=ax)
# ax.set_title('Cost Breakdown', fontsize=24)
# ax.set_ylabel('Costs in kCHF/year', fontsize=24)
# # ax.set_xlabel('Scenario', fontsize=24)
# plt.xticks(rotation=0, fontsize=24)  # Rotate x-axis labels
# plt.yticks(fontsize=24)

# # Positioning the legend on the right side of the plot
# ax.legend(title='Cost Category', loc='upper left', bbox_to_anchor=(1, 1), title_fontsize='24', fontsize=24)

# plt.tight_layout()
# plt.savefig('Grid-Charges.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# IMPACT OF WASTE HEAT RECOVERY
#------------------------------------------------------------------------------
# # Data based on your description, structured as a dictionary
# data = {
#     'Installation': [202.67, 193.32],
#     'Operation': [66.05, 88.61],
#     'Maintenance': [11.69, 10.88],
#     'Revenues WHR': [-16.47, 0.00]  # Negative values as required
# }

# # Creating DataFrame
# df = pd.DataFrame(data, index=['SS75', 'w/o WHR'])

# # Setting seaborn style and color palette
# sns.set(style="white", palette="husl")
# sns.set_context("talk", font_scale=1.2)  # Adjusting font scale globally

# # Plotting the updated stacked bar chart
# fig, ax = plt.subplots(figsize=(10, 10))

# # Using DataFrame plot method to handle stacking directly
# df.plot(kind='bar', stacked=True, ax=ax, color=sns.color_palette("husl"))

# # Customizing plot features
# ax.set_title('Cost Breakdown', fontsize=24)
# ax.set_ylabel('Costs in kCHF/year', fontsize=24)
# # ax.set_xlabel('Scenario', fontsize=24)
# plt.xticks(rotation=0, fontsize=24)  # Set x-axis labels to horizontal
# plt.yticks(fontsize=24)

# # Positioning the legend on the right side of the plot
# ax.legend(title='Cost Category', loc='upper left', bbox_to_anchor=(1, 1), title_fontsize='24', fontsize=24)

# plt.tight_layout()
# plt.savefig('impactWHR.png', bbox_inches='tight', dpi=150)
# plt.show()

#------------------------------------------------------------------------------
# PWA
#------------------------------------------------------------------------------

# # Path to the Excel file
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-12_23-02-37 - ultimate - SS90 - with PWA.xlsx'

# # Load the data from the specified columns in the Excel file
# df_optimization_results = pd.read_excel(file_path, sheet_name='PWA', usecols=['P_ELY', 'P_ELY_PWA', 'ETA_ELY', 'P_FC_in', 'P_FC', 'ETA_FC'])

# # Constants provided by the user

# S_ELY = 330726.362
# S_FC = 226129.5564


# # Calculate the ratios for x and y axis
# df_optimization_results['Input Power'] = df_optimization_results['P_ELY'] / S_ELY
# df_optimization_results['Output Power'] = df_optimization_results['P_ELY_PWA'] / S_ELY

# # Set up the plot with increased font size
# fig, ax = plt.subplots(figsize=(8, 8))
# sc = ax.scatter(df_optimization_results['Input Power'], df_optimization_results['Output Power'], 
#                 c=df_optimization_results['ETA_ELY'], cmap='viridis', s=250)

# # Color bar setup
# cbar = plt.colorbar(sc)
# cbar.set_label('Efficiency in %', fontsize=22)

# # Labeling the plot with increased font size
# ax.set_xlabel('Input Power', fontsize=22)
# ax.set_ylabel('Output Power', fontsize=22)
# ax.set_title('Electrolyser Input/Output Relation Normalised', fontsize=24)

# # Increase the font size for the tick labels
# ax.tick_params(axis='both', which='major', labelsize=20)
# cbar.ax.tick_params(labelsize=20)

# # Display the plot
# plt.tight_layout()
# # plt.savefig('ELY-PWA.png', bbox_inches='tight', dpi=150)
# plt.show()

# #------------------------------------------------------------------------------

# # Calculate the ratios for x and y axis using new columns
# df_optimization_results['Input Power FC'] = df_optimization_results['P_FC_in'] / S_FC
# df_optimization_results['Output Power FC'] = df_optimization_results['P_FC'] / S_FC

# # Set up the plot for Fuel Cell data
# fig, ax = plt.subplots(figsize=(8, 8))
# sc = ax.scatter(df_optimization_results['Input Power FC'], df_optimization_results['Output Power FC'], 
#                 c=df_optimization_results['ETA_FC'], cmap='viridis', s=250)

# # Color bar setup
# cbar = plt.colorbar(sc)
# cbar.set_label('Efficiency in %', fontsize=22)

# # Labeling the plot with increased font size
# ax.set_xlabel('Input Power', fontsize=22)
# ax.set_ylabel('Output Power', fontsize=22)
# ax.set_title('Fuel Cell Input/Output Relation Normalised', fontsize=24)

# # Increase the font size for the tick labels
# ax.tick_params(axis='both', which='major', labelsize=20)
# cbar.ax.tick_params(labelsize=20)

# # Display the plot
# plt.tight_layout()
# #plt.savefig('FC-PWA.png', bbox_inches='tight', dpi=150)
# plt.show()


#------------------------------------------------------------------------------
# HEAT
#------------------------------------------------------------------------------

# # Path to the uploaded Excel file
# file_path = 'C:/Users/fism/Desktop/MA_thesis/02_modeling_and_optimization/results/optimization_results_2024-05-07_09-58-27 - current - SS75.xlsx'

# # Load the data from the Excel file
# df_optimization_results = pd.read_excel(file_path, sheet_name='Optimization Results')

# # Assuming the dataset is hourly and covers a full year (8760 hours)
# # We need to create a datetime index for the year 2019 (or the relevant year)
# date_range = pd.date_range(start='2019-01-01', periods=8760, freq='H')
# df_optimization_results['timestamp'] = date_range
# df_optimization_results.set_index('timestamp', inplace=True)

# # Convert power values from Watts to kilowatts

# df_optimization_results[['P_th_HT', 'P_th_MT',
#     'z1_35degC_kWh', 'z1_60degC_kWh',
#     'z2_35degC_kWh', 'z2_60degC_kWh',
#     'z3_35degC_kWh', 'z3_60degC_kWh']] /= 1000

# # Resample data to daily sums

# daily_data = df_optimization_results[['P_th_HT', 'P_th_MT',
#     'z1_35degC_kWh', 'z1_60degC_kWh',
#     'z2_35degC_kWh', 'z2_60degC_kWh',
#     'z3_35degC_kWh', 'z3_60degC_kWh']].resample('D').sum()

# # Set the seaborn style and palette
# sns.set(style="white")
# colors = sns.color_palette("Set2")

# # Increase font sizes globally
# plt.rcParams.update({'font.size': 20, 'axes.titlesize': 22, 'axes.labelsize': 22, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'legend.fontsize': 18})

# # Assuming 'daily_data' is a DataFrame already defined and containing the necessary data
# fig, ax = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

# # Calculate sums and convert to MWh or MW
# P_th_HT_sum = daily_data['P_th_HT'].sum() / 1000  # Convert to MW
# z3_60degC_kWh_sum = daily_data['z3_60degC_kWh'].sum() / 1000  # Convert to MWh
# P_th_MT_sum = daily_data['P_th_MT'].sum() / 1000  # Convert to MW
# z3_35degC_kWh_sum = daily_data['z3_35degC_kWh'].sum() / 1000  # Convert to MWh

# # First plot with specified colors and labels including sums
# ax[0].fill_between(daily_data.index, daily_data['P_th_HT'], label=f'Recovered Heat, Total: {P_th_HT_sum:.2f} MWh', color=colors[0], alpha=0.8)
# ax[0].plot(daily_data.index, daily_data['z3_60degC_kWh'], label=f'Zone 3 Demand, Total: {z3_60degC_kWh_sum:.2f} MWh', color=colors[1])
# ax[0].set_title('Domestic Hot Water')
# # ax[0].set_xlabel('Date')
# ax[0].set_ylabel('Power in kW')
# ax[0].legend()

# # Second plot with specified colors and labels including sums
# ax[1].fill_between(daily_data.index, daily_data['P_th_MT'], label=f'Recovered Heat, Total: {P_th_MT_sum:.2f} MWh', color=colors[2], alpha=0.8)
# ax[1].plot(daily_data.index, daily_data['z3_35degC_kWh'], label=f'Zone 3 Demand, Total: {z3_35degC_kWh_sum:.2f} MWh', color=colors[3])
# ax[1].set_title('Space Heating')
# ax[1].set_xlabel('Date')
# ax[1].set_ylabel('Power in kW')
# ax[1].legend(loc='upper center')

# plt.tight_layout()
# plt.savefig('WHR - Zone 3.png', bbox_inches='tight', dpi=150)
# plt.show()
