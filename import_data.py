# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:01:35 2024
@author: Maxime Fischer 

This module imports the solar irradiance for Rubigen from the 
Rubigen_2019-2022_irradiance_Hourly.xlsx file as well as the power demand from 
the tesla supercharger station from 200923_Lastang 2015-2020_Rubigen.xlsx
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io

# Set seaborn color palettes
sns.set_palette("Accent")
palette = sns.color_palette()


# input_path  = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx'
# demand_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\demand_data\200923_Lastang 2015-2020_Rubigen.xlsx'

class WeeklyTimeline:
    def __init__(self, input_path, demand_path, season):
        self.input_path = input_path
        self.demand_path = demand_path
        self.season = season  # 'summer' or 'winter'
        self.irradiance = None
        self.T_amb = None
        self.P_demand = None
        self.index = None

    def import_data(self):
        sheet_name = '2019_week_sum' if self.season == 'summer' else '2019_week_win'
        df_input = pd.read_excel(self.input_path, sheet_name=sheet_name)
        self.irradiance = df_input['ALLSKY_SFC_SW_DWN [Wh/m2]'].values   
        self.T_amb = df_input['T2M [°C]'].values
        
        df_demand = pd.read_excel(self.demand_path, sheet_name=sheet_name)
        df_demand = pd.read_excel(self.demand_path, sheet_name='30%BEV+heavyFri&Sat1WeekOnly')

        self.P_demand = df_demand.groupby(df_demand.index // 4)['Wirkleistung [kW]'].mean().values
        self.P_demand = (self.P_demand * 1000)  # Power demand in [W]
        
        self.index = range(len(self.irradiance))
        
        # Add logic to import heating demand data
        heat_sheet_name = 'heat_week_sum' if self.season == 'summer' else 'heat_week_win'
        df_heat_demand = pd.read_excel(self.demand_path, sheet_name=heat_sheet_name)
        
        return self.irradiance, self.P_demand, self.T_amb, df_input, df_demand, df_heat_demand
    
    
class MonthlyTimeline:
    def __init__(self, input_path, demand_path, season):
        self.input_path = input_path
        self.demand_path = demand_path
        self.season = season  # 'summer' or 'winter'
        self.irradiance = None
        self.T_amb = None
        self.P_demand = None
        self.index = None

    def import_data(self):
        sheet_name = '2019_month_sum' if self.season == 'summer' else '2019_month_win'
        df_input = pd.read_excel(self.input_path, sheet_name=sheet_name)
        self.irradiance = df_input['ALLSKY_SFC_SW_DWN [Wh/m2]'].values
        self.T_amb = df_input['T2M [°C]'].values
        
        df_demand = pd.read_excel(self.demand_path, sheet_name=sheet_name)
        self.P_demand = df_demand.groupby(df_demand.index // 4)['Wirkleistung [kW]'].mean().values
        self.P_demand = (self.P_demand * 1000)  # Power demand in [W]
        
        self.index = range(len(self.irradiance))
    
        heat_sheet_name = 'heat_month_sum' if self.season == 'summer' else 'heat_month_win'
        df_heat_demand = pd.read_excel(self.demand_path, sheet_name=heat_sheet_name)
        
        return self.irradiance, self.P_demand, self.T_amb, df_input, df_demand, df_heat_demand


class YearlyTimeline:
    def __init__(self, input_path, demand_path, year, demand_choice='Normal'):
        self.input_path = input_path
        self.demand_path = demand_path
        self.year = year  # Added to specify the year
        self.demand_choice = demand_choice
        self.irradiance = None
        self.T_amb = None
        self.P_demand = None

    def import_data(self):
        # Read the specified year's sheet and the specific column for irradiance
        df_input = pd.read_excel(self.input_path, sheet_name=str(self.year))
        self.irradiance = df_input['ALLSKY_SFC_SW_DWN [Wh/m2]'].values  # yearly solar irradiance in [Wh/m2]
        self.T_amb = df_input['T2M [°C]'].values
        
        # Choose the appropriate sheet for demand data based on user's choice if year is 2019
        if self.year == '2019':
            if self.demand_choice == 'Normal':
                df_demand = pd.read_excel(self.demand_path, sheet_name='2019')
            elif self.demand_choice == 'Increased':
                df_demand = pd.read_excel(self.demand_path, sheet_name='30%BEV+heavyFri&Sat')
        else:
            df_demand = pd.read_excel(self.demand_path, sheet_name='2019')  # Default to '2019' sheet for other years
        
        # Aggregate demand from 15min to 1h timesteps by means
        self.P_demand = df_demand.groupby(df_demand.index // 4)['Wirkleistung [kW]'].mean().values
        self.P_demand = (self.P_demand * 1000)  # Power demand in [W]

        # Create an index based on the length of the data
        self.index = range(len(self.irradiance))

        df_heat_demand = pd.read_excel(self.demand_path, sheet_name='heat_year')
        
        return self.irradiance, self.P_demand, self.T_amb, df_input, df_demand, df_heat_demand

# Declare as global variables
irradiance = None
P_demand   = None


def get_data(input_path, demand_path):
    global irradiance, P_demand
    timeline_choice = input("Enter 'week', 'month' or 'year' for the timeline: ")
    season_choice = None  # Initialize season_choice to ensure it's always defined
    
    if timeline_choice.lower() == 'week' or timeline_choice.lower() == 'month':
        season_choice = input("Do you want to consider a 'summer' or 'winter' season? ")
        while season_choice not in ['summer', 'winter']:
            print("Invalid season choice. Please choose 'summer' or 'winter'.")
            season_choice = input("Do you want to consider a 'summer' or 'winter' season? ")
        
        if timeline_choice.lower() == 'week':
            weekly_timeline = WeeklyTimeline(input_path, demand_path, season_choice)
            irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand = weekly_timeline.import_data()
            plot_data(weekly_timeline, f'{season_choice.capitalize()} Week 2019')
            
        elif timeline_choice.lower() == 'month':
            monthly_timeline = MonthlyTimeline(input_path, demand_path, season_choice)
            irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand = monthly_timeline.import_data()
            plot_data(monthly_timeline, f'{season_choice.capitalize()} Month 2019')
    
    elif timeline_choice.lower() == 'year':
        valid_years = ['2019', '2020', '2021', '2022']
        year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")
        while year_choice not in valid_years:
            print("Invalid year choice. Please choose from 2019, 2020, 2021, or 2022.")
            year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")

        

        demand_choice = 'Normal'  # Default to normal unless specified otherwise
        if year_choice == '2019':
            demand_choice = input("Enter 'Normal' or 'Increased' for the demand scenario: ")
            while demand_choice not in ['Normal', 'Increased']:
                print("Invalid demand choice. Please choose 'Normal' or 'Increased'.")
                demand_choice = input("Enter 'Normal' or 'Increased' for the demand scenario: ")


        yearly_timeline = YearlyTimeline(input_path, demand_path, year_choice, demand_choice)
        irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand = yearly_timeline.import_data()
        
    # elif timeline_choice.lower() == 'year':
    #     valid_years = ['2019', '2020', '2021', '2022']
    #     year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")
    #     while year_choice not in valid_years:
    #         print("Invalid year choice. Please choose from 2019, 2020, 2021, or 2022.")
    #         year_choice = input("Enter the year you want to consider (2019, 2020, 2021, 2022): ")
            
    #     yearly_timeline = YearlyTimeline(input_path, demand_path, year_choice)
    #     irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand = yearly_timeline.import_data()
       
        plot_data(yearly_timeline, f'{year_choice}')
        
    else:
        print("Invalid timeline choice. Please enter 'week', 'year', or 'month'.")
    
    return irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand, timeline_choice, season_choice

    

if __name__ == "__main__":
    get_data()

#------------------------------------------------------------------------------
# PLOTTING
#------------------------------------------------------------------------------

# DON'T DELETE THIS ALTERNATIVE PLOT!!!!---------------------------------------
# def plot_data(timeline, title_prefix):
#     # Create a figure with two subplots (one above the other)
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 8))

#     # Define a common font size
#     common_fontsize = 24  # Adjust the size as needed

#     # Plot Power Demand on ax1
#     ax1.plot(timeline.index, timeline.P_demand / 1000, label='Power Demand', color='blue')
#     ax1.set_ylabel('kW', fontsize=common_fontsize)
#     ax1.legend(loc='upper left', fontsize=common_fontsize)
#     ax1.tick_params(axis='both', which='major', labelsize=common_fontsize)

#     # Plot Solar Irradiance on ax2
#     ax2.plot(timeline.index, timeline.irradiance / 1000, label='Irradiance', color='orange')
#     ax2.set_xlabel('Time (hours)', fontsize=common_fontsize)
#     ax2.set_ylabel('kW/m2', fontsize=common_fontsize)
#     ax2.legend(loc='upper right', fontsize=common_fontsize)
#     ax2.tick_params(axis='both', which='major', labelsize=common_fontsize)

#     # Set the figure title
#     fig.suptitle(f'{title_prefix}: Input data', fontsize=common_fontsize + 4)  # Slightly larger font size for title

#     plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the layout to make room for the figure title
#     plt.show()

# DON'T DELETE THIS ALTERNATIVE PLOT!!!!!!-------------------------------------
def plot_data(timeline, title_prefix):
    fig, ax1 = plt.subplots(figsize=(20, 6))

    # Plot Power Demand on the left y-axis
    ax1.fill_between(timeline.index, timeline.P_demand / 1000, label='Power Demand', color=palette[4], alpha=0.9)
    ax1.set_xlabel('Time (hours)', fontsize=24)
    ax1.set_ylabel('Power Demand (kW)', fontsize=24)
    #ax1.legend(loc='upper left', fontsize=24)
    ax1.tick_params(axis='both', which='major', labelsize=24)
    
    ax1.legend(loc='lower left', bbox_to_anchor=(0, -0.3),
          fancybox=False, shadow=False, ncol=1, fontsize=24)

    # Create a second y-axis for Solar Irradiance on the right side
    ax2 = ax1.twinx()
    ax2.fill(timeline.index, timeline.irradiance / 1000, label='Irradiance', color=palette[3], alpha=0.7)
    ax2.set_ylabel('Irradiance (kW/m2)', fontsize=24)
    #ax2.legend(loc='upper right', fontsize=24)
    ax2.tick_params(axis='both', which='major', labelsize=24)
    
    ax2.legend(loc='lower right', bbox_to_anchor=(1, -0.3),
          fancybox=False, shadow=False, ncol=1, fontsize=24)

    # Set titles and labels
    fig.suptitle(f'{title_prefix} Power Demand and Solar Irradiance', fontsize=28)
    plt.grid(False)
    plt.tight_layout()
    plt.show()