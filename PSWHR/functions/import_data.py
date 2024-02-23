# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:01:35 2024

@author: fism
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set seaborn color palettes
sns.set_palette("Set2")
palette = sns.color_palette()

class WeeklyTimeline:
    def __init__(self, irradiance_path, demand_path):
        self.irradiance_path = irradiance_path
        self.demand_path = demand_path
        self.irradiance = None
        self.P_demand = None

    def import_data(self):
        df_input = pd.read_excel(self.irradiance_path)   # read irradiance file
        self.irradiance = df_input['G'].values           # hourly solar irradiance in [W/m2]
        df_demand = pd.read_excel(self.demand_path)      # read demand file

        # Aggregate demand from 15min to 1h timesteps by summing every 4 values
        self.P_demand = df_demand.groupby(df_demand.index // 4)['power_demand_kW'].mean().values
        self.P_demand = (self.P_demand * 1000)  # Power demand in [W] 

        # Create an index based on the length of the data
        self.index = range(len(self.irradiance))

        return self.irradiance, self.P_demand, df_input, df_demand

# class YearlyTimeline:
#     def __init__(self, irradiance_path, demand_path):
#         self.irradiance_path = irradiance_path
#         self.demand_path = demand_path
#         self.irradiance = None
#         self.P_demand = None

#     def import_data(self):
#         # Read the '2019' sheet and the specific column for irradiance
#         df_input = pd.read_excel(self.irradiance_path, sheet_name='2019')
#         self.irradiance = df_input['ALLSKY_SFC_SW_DWN [Wh/m2]'].values  # yearly solar irradiance in [Wh/m2]
#         df_demand = pd.read_excel(self.demand_path)

#         # Aggregate demand from 15min to 1h timesteps by means
#         self.P_demand = df_demand.groupby(df_demand.index // 4)['Wirkleistung [kW]'].mean().values
#         self.P_demand = (self.P_demand * 1000)  # Power demand in [W]

#         # Create an index based on the length of the data
#         self.index = range(len(self.irradiance))

        return self.irradiance, self.P_demand, df_input, df_demand

class YearlyTimeline:
    def __init__(self, irradiance_path, demand_path, year):
        self.irradiance_path = irradiance_path
        self.demand_path = demand_path
        self.year = year  # Added to specify the year
        self.irradiance = None
        self.P_demand = None

    def import_data(self):
        # Read the specified year's sheet and the specific column for irradiance
        df_input = pd.read_excel(self.irradiance_path, sheet_name=str(self.year))
        self.irradiance = df_input['ALLSKY_SFC_SW_DWN [Wh/m2]'].values  # yearly solar irradiance in [Wh/m2]
        df_demand = pd.read_excel(self.demand_path)

        # Aggregate demand from 15min to 1h timesteps by means
        self.P_demand = df_demand.groupby(df_demand.index // 4)['Wirkleistung [kW]'].mean().values
        self.P_demand = (self.P_demand * 1000)  # Power demand in [W]

        # Create an index based on the length of the data
        self.index = range(len(self.irradiance))

        return self.irradiance, self.P_demand, df_input, df_demand


def plot_data(timeline, title_prefix):
    fig, ax1 = plt.subplots(figsize=(15, 8))

    # Plot Power Demand on the left y-axis
    ax1.plot(timeline.index, timeline.P_demand / 1000, label='Power Demand (kW)', color=palette[0])
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Power Demand (kW)')
    ax1.legend(loc='upper left')

    # Create a second y-axis for Solar Irradiance on the right side
    ax2 = ax1.twinx()
    ax2.plot(timeline.index, timeline.irradiance / 1000, label='Irradiance (kW/m2)', color=palette[1])
    ax2.set_ylabel('Irradiance (kW/m2)')
    ax2.legend(loc='upper right')

    # Fill the area between P_demand and irradiance where P_demand is above irradiance
    mask = timeline.P_demand > timeline.irradiance
    ax1.fill_between(timeline.index, timeline.P_demand / 1000, timeline.irradiance / 1000, where=mask,
                     color=palette[3], alpha=0.5)

    # Set titles and labels
    fig.suptitle(f'{title_prefix} Power Demand and Solar Irradiance')
    plt.tight_layout()
    plt.show()


# Declare as global variables
irradiance = None
P_demand   = None

# def get_data():
#     global irradiance, P_demand
#     timeline_choice = input("Enter 'week' or 'year' for the timeline: ")
#     if timeline_choice.lower() == 'week':
#         weekly_timeline = WeeklyTimeline(
#             r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\input_data\Inputs_BoundaryLoads_week.xlsx',
#             r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\demand_data\demand_data_weekly.xlsx'
#         )
        
#         # Save variables in the main script
#         df_input, df_demand, irradiance, P_demand = weekly_timeline.import_data()
#         plot_data(weekly_timeline, 'Weekly')
        
#     # elif timeline_choice.lower() == 'year':
#     #     yearly_timeline = YearlyTimeline(
#     #         r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\input_data\Inputs_BoundaryLoads_year_2022.xlsx',
#     #         r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\demand_data\demand_data_yearly.xlsx'
#     #     )
#     elif timeline_choice.lower() == 'year':
#         yearly_timeline = YearlyTimeline(
#         r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
#         r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\demand_data\demand_data_yearly.xlsx'
#         )
        
#         # Save variables in the main script
#         df_input, df_demand, irradiance, P_demand = yearly_timeline.import_data()
#         plot_data(yearly_timeline, 'Yearly')
#     else:
#         print("Invalid timeline choice. Please enter 'week' or 'year'.")
#         exit()

def get_data():
    global irradiance, P_demand
    timeline_choice = input("Enter 'week' or 'year' for the timeline: ")
    if timeline_choice.lower() == 'week':
        weekly_timeline = WeeklyTimeline(
            # r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\input_data\Inputs_BoundaryLoads_week.xlsx',
            # r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\demand_data\demand_data_weekly.xlsx'
            r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\Inputs_BoundaryLoads_week.xlsx',
            r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\demand_data_weekly.xlsx'
        )
        
        # Save variables in the main script
        df_input, df_demand, irradiance, P_demand = weekly_timeline.import_data()
        plot_data(weekly_timeline, 'Weekly')
        
        
    elif timeline_choice.lower() == 'year':
        # Ask the user for the year they wish to consider
        valid_years = ['2019', '2020', '2021', '2022']
        year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")
        while year_choice not in valid_years:
            print("Invalid year choice. Please choose from 2019, 2020, 2021, or 2022.")
            year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")
        
        yearly_timeline = YearlyTimeline(
            r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
            r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\demand_data_yearly.xlsx',
            year_choice  # Pass the chosen year to the YearlyTimeline class
        )
        
        # Save variables in the main script
        df_input, df_demand, irradiance, P_demand = yearly_timeline.import_data()
        plot_data(yearly_timeline, f'Yearly - {year_choice}')
    else:
        print("Invalid timeline choice. Please enter 'week' or 'year'.")
        exit()

    # Now you can use irradiance and P_demand in the main script
    print("Irradiance:", irradiance)
    print("Power Demand:", P_demand)
    
    return irradiance, P_demand, df_input, df_demand

if __name__ == "__get_data__":
    get_data()