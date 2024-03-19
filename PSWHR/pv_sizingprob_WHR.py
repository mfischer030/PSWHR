# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:26:35 2023

@author: fism
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import Model, GRB
import sys # 
from dash import Dash, html, dcc

"""
Introduction: sizing and operational optimization problem
This model solves a sizing and operational problem with a MILP approach for 
a Hydrogen energy storage system sonsisting of
- An electrolyser (ELY)
- A compressor (C)
- A pressurized hydrogen storage tank (TANK)
- A fuel cell (FC)

The HESS is feeded by an PV-field (PV) and is connected (or not) to the grid.
The objective of this optimization is to minimze the cost function:
    cost = cost_inst + cost_op + cost_maint + cost_start
"""
#------------------------------------------------------------------------------
# Choose user or adapt paths in config.py module
#------------------------------------------------------------------------------
user = 'maxime'    # 'christian', 'gabriele', 'maxime'

from config import paths_configuration
input_path, demand_path, heat_path, function_path, export_path = paths_configuration(user)
#------------------------------------------------------------------------------
# Setting up the model
#------------------------------------------------------------------------------
# Choose energy tariff
energy_tariff = "Blue" # Choose from "Green","Blue","Grey"

# Flag to include/exclude battery in the optimization model
include_battery = False

#------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------
# path to the functions directory => change in config.py if needed
sys.path.append(function_path) 

# Import data and generate variables using the get_data function---------------
from import_data import get_data
irradiance, P_demand, T_amb, df_input, df_demand, timeline_choice = get_data(input_path, demand_path)

# Heat demand data (from NEST building @ empa, year 2022)----------------------
df_demand_heat = pd.read_excel(heat_path)

heat_35degC_demand = df_demand_heat['Heating_35degC_kW'].values  
heat_65degC_demand = df_demand_heat['DHW_65degC_kW'].values     

from plotting_module import heat_demand_plot
heat_demand_plot(heat_35degC_demand,heat_65degC_demand)
#------------------------------------------------------------------------------
# Input parameters
#------------------------------------------------------------------------------
# General parameters
k    = 1.4                               # Ratio cp/cv [-]
R_H2 = 4.1242 * 1000                     # Individual Gas constant H2 [J/kg*K]
HHV  = 39.39 * 3600 * 1000               # Higher heating value of H2 in [J/kg] = 39.39 kWh/kg

project_lifetime     = 25   # in [y] 2023_Wang et al - 20 years
annual_interest_rate = 0.04 # Discount rate, as encouraged by EU, from Rox

# Electricity prices in [EUR/MWh]
cost_imp_el = df_input['price_Eur_MWh'].values # hourly cost to import 
cost_exp_el = df_input['Price_DayAhed'].values # hourly cost to export 

# Revenues from WHR in [€/kWh]
cost_export_heatMT = 0.2 #0.09 
cost_export_heatHT = 0.2 #0.1189

# Time parameters
nHours = len(df_input)                        # number of hours simulated
Time   = np.arange(1, nHours + 1)             # time vector
days   = nHours / 24                          # number of days
months = df_input['MO'].unique()              # Unique months in the DataFrame
weeks  = days / 7                             # number of weeks
deltat = 3600                                 # time step (s)

# Convert nHours to datetime format
if timeline_choice == 'year':
    df_input['ts'] = pd.Timestamp('2019-01-01') + pd.to_timedelta(Time - 1, unit='H')
    print(df_input['ts'].iloc[0])  # Print first datetime to verify

#------------------------------------------------------------------------------
# Efficiencies of components in [-]
#------------------------------------------------------------------------------
# PV: from Roxanne
# ELY: 2023_van_der_roest (74-79%)
# C: from Roxanne | 2020_Pan et al: 0.9
# TANK: 2023_Wang et al
# FC: 2021_cigolotti Comprehensive Review on FC Technology for Stationary Applications: Electric Efficiency PEMFC: [38,38,37,40] in [%]

eta = {'PV': 0.21,'ELY': 0.74,'C': 0.7526,'TANK': 0.95,'FC': 0.5,}              # Realistic
#eta = {'PV': 0.21,'ELY': 0.87,'C': 0.8763,'TANK': 0.975,'FC': 0.75}
#eta = {'PV': 0.21,'ELY': 1,'C': 1,'TANK': 1,'FC': 1}                           # Optimal      

#------------------------------------------------------------------------------
# Unit prices of components / capital costs in [€/W] 
#------------------------------------------------------------------------------
# PV:   in [€/W] | 2023_Tay Son Le: 881 USD/year | 2018_Gabrielli 300€/m2 
# BAT:  Unit price of the battery in CHF/kWh
# ELY:  in [€/W] | 2023_IEA-GlobalH2Review: PEM - 2kUSD/kW - reduction to 600 USD/kW
# C:    in [€/W] | 2020_Pan et al: 1228 (¥/kW)
# TANK: in [EUR/J], 1644 [EUR/kgH2] | 2018_Gabrielli: {20.7;13.6;10.9} in [€/kWh] | 2023_A Review on the Cost Analysis of H2 Gas Storage Tanks for FCV: 2020:9.34 €/kWh 2025: 8.40€/kWh Ultimate: 7.47€/kWh
# FC:   in [€/W] (Wang et. al 2023 - 2kUSD/kW) | 2018_Gabrielli {2160;1680;1320} [€/kW]
# HEX: in [€/m2]
Fixed_HEX = 5291.9 # Fixed cost for heat exchanger [EUR]

# UP = {'PV': 0.8,'BAT': 500,'ELY': 1.7,'C': 0.0076,'TANK': 9.34/(3.6*10**6),'FC': 2, 'HP': 0.576, 'HEX': 77.79}       # Realistic
#UP = {'PV': 0.8,'BAT': 500,'ELY': 1.128,'C': 0.00506,'TANK': 8.4/(3.6*10**6),'FC': 1.3, 'HP': 0.238, 'HEX': 40}
UP  = {'PV': 0.8, 'BAT': 500,'ELY': 0.556,'C': 0.0038,'TANK': 7.47/(3600000),'FC': 0.6, 'HP': 0, 'HEX': 0}         # Optimal


# UP = {'PV': 0.8,'ELY': 1.7/10,'C': 0.0076,'TANK': 1644/(10*HHV),'FC': 1.680/10} #Old values

# Annual maintenance cost as fraction of total cost => from Roxanne------------
maintenance = {
    'PV': 0.0158,       # Annual maintenance PV
    'BAT': 0.02,        # Annual cost maintenance battery
    'ELY': 0.02,        # Annual maintenance electrolyser
    'C': 0.08,          # Annual maintenance compressor
    'TANK': 0.03,       # Annual maintenance storage tank
    'FC': 0.02,         # Annual maintenance fuel cell
    'HP': 0.015,        # Annual maintenance heat pump
    'HEX': 0.01         # Annual cost maintenance HEX, frac total ann cost
}

# Lifetime of the components in [years]----------------------------------------
life = {'PV': 25,    # Roxanne: 30 | 2023_Tya Son Le: 25 years
        'BAT': 10,   # From Gabriele
        'ELY': 10,   # 2023_Wang et al: 5 years | 2023_Tay Son Le: 15 years || Maxime: Before value was 20
        'C': 20,     # 2020_Pan et al: 20 years
        'TANK': 35,  # 2023_Wang et al: 20 years | 2023_Tay Son Le: 25 years
        'FC': 5,     # 2023_Wang et al: 5 years | 2023_Tay Son Le: 5 years || Maxime: Before value was 20
        'HP': 20,
        'HEX': 20    # from Roxanne
}

# Daily energy demand in [Wh] - P_demand is in [W]
E_demand_day = sum(P_demand)/days # Gabriele: this is [Wh], correct, Maxime: yes

# Maximal daily power demand in [W]
P_peak_max = np.max(P_demand)

# Battery Parameters
battery_params = {
    'eff_ch': 0.95,         # Battery charging efficiency
    'eff_disch': 0.95,      # Battery discharging efficiency
    'eff_sd': 0.99,         # Battery self-discharge efficiency
    'C_b_max': 50*3600000,  # Max battery capacity in [J]
    'C_b_min': 0,           # Min battery capacity in [J] - typically set to 0
    'SOC_max': 0.8,         # Max state of charge of the battery to increase lifetime
    'SOC_min': 0.2          # Min state of charge of the battery to increase lifetime
}

#------------------------------------------------------------------------------
# Defining maximal SIZES of the components
#------------------------------------------------------------------------------

# Area_PV_max  = P_peak_max / (eta['PV'] * np.mean(irradiance)) # Maximum PV area [m2]
Area_PV_max  = 3500 # Maximum PV area [m2] => how to set, quantify the maximum pv area?

# Electrolyser max and min nominal power (W)   
S_ELY_max = 1000*1000   # Maximal Electrolyzer size in [W]
S_ELY_min = 0           # Min. size ELY where problem is feasible [W] - from Rox

# Calculating the spezific work of the compresso, from Minutillo et al. 2021
# (Analyzing the levelized cost of hydrogen eq 1+2) => from Roxanne 
T_in_H2 = 65 + 273.15                       # H2 temperature (=T_cat=T_an) [K]  
p_out   = 350                               # Compressor outlet pressure [bar] = H2 storage pressure 
p_in    = 30                                # Compressor inlet pressure [bar]  PEM electrolyzer at Empa works at 30 bar 
L_is_C  = (k/(k-1)) * R_H2 * T_in_H2 * (((p_out/p_in)**((k-1)/k)) - 1) # Specific work compressor [J/kg] 

# Maximal TANK energy capacity (J) => 4 days storage capacity
S_TANK_max = E_demand_day * 14 * 3600                     # E_demand_day in [Wh]  
S_TANK_H2_max = S_TANK_max / HHV                          # equivalent in kg_H2

# Maximal FC size as the maximal power demand divided by eff in [W]
#S_FC_max = P_peak_max / eta["FC"]
S_FC_max = 1000*1000  
"""
Range PEMFC = [10W;1MW] from 2021_cigolotti Comprehensive Review on Fuel Cell 
Technology for Stationary Applications
"""
# Waste Heat Recovery ---------------------------------------------------------

T_in  = 57    # Inlet temperature cooling water to HEX in [°C]
T_out = 62    # Temperature cooling water to applications [°C]
T_HEX = 64    # Outlet temperature cooling water from HEX [°C]
c_p   = 4186  # Specific heat capacity cooling water [J/kg/K]

P_th_ELY_max = (1 - eta['ELY']) * S_ELY_max # max heat recovered from ELY in [W]
P_th_FC_max  = (1 - eta['FC'])  * S_FC_max  # max heat recovered from FC in [W]
P_th_max = P_th_ELY_max + P_th_FC_max

m_cw_FC_max  = P_th_FC_max / (c_p * (T_out - T_in))      # Cooling water maximum mass flow in [kg/s]
m_cw_ELY_max = P_th_ELY_max / (c_p * (T_out - T_in))     # Cooling water maximum mass flow in [kg/s]
m_cw_max     = P_th_max / (c_p * (T_out - T_in))         # Cooling water maximum mass flow in [kg/s]

T_MT_in   = 26        # MT DHN water inlet temperature HEX [°C] from NEST
T_MT_out  = 36        # MT DHN water outlet temperature HEX [°C] from NEST
T_HT_out  = 66        # HT DHN water outlet temperature HP [°C] from NEST (domestic hot water)
T_log     = ((T_out - T_MT_out) - (T_in - T_MT_in)) / np.log((T_out - T_MT_out) / (T_in - T_MT_in)) # Logarithmic mean temp difference in HEX
U_HEX     = 2         # Overall heat transfer coeff HEX [kW/m2/K], swedish thesis

S_HEX_ELY_max = P_th_ELY_max / (U_HEX * T_log)  # Maximum heat exchanger surface [m2]
S_HEX_FC_max  = P_th_FC_max / (U_HEX * T_log)   # Maximum heat exchanger surface [m2]
S_HEX_max     = S_HEX_ELY_max + S_HEX_FC_max

# Coefficient of performance (COP)
COP_carnot= T_HT_out / (T_HT_out - T_out)             # Maximum COP HP
COP       = 0.5 * COP_carnot                          # Real COP, as in Tiktak
#------------------------------------------------------------------------------
# Prompt the user to choose "Grid" or "Off-Grid" scenario
scenario_choice = input("Enter 'grid' for grid-connected scenario or 'off-grid' for an off-grid scenario: ").lower()
# Check the user's choice and set relevant parameters accordingly
if scenario_choice == 'grid':
    P_imp_ub = 1000000                 # Upper bound for P_imp in grid-connected scenario, Trafo limit is 1MW
    P_exp_ub = 1000000                 # Upper bound for P_exp in grid-connected scenario, Trafo limit is 1MW
elif scenario_choice == 'off-grid':
    P_imp_ub = 0                       # No imported power in off-grid scenario
    P_exp_ub = 0                       # No exported power in off-grid scenario
else:
    print("Invalid scenario choice. Please enter 'grid' or 'off-grid'.")
    sys.exit(1)  # Exit the program in case of an invalid choice
    
#------------------------------------------------------------------------------
# Define the optimization problem
m = Model()
#------------------------------------------------------------------------------
# Design variables-------------------------------------------------------------
# 1D-Variables
Area_PV = m.addVar(lb=0,         ub=Area_PV_max, name='Area_PV')
S_ELY   = m.addVar(lb=S_ELY_min, ub=S_ELY_max,   name='S_ELY')
S_TANK  = m.addVar(lb=0,         ub=S_TANK_max,  name='S_TANK')
S_FC    = m.addVar(lb=0,         ub=S_FC_max,    name='S_FC')

# Time dependent variables
P_ELY   = m.addVars(nHours, lb=0, ub=S_ELY_max,  name='P_ELY')
E_TANK  = m.addVars(nHours, lb=0, ub=S_TANK_max, name='E_TANK')
P_FC    = m.addVars(nHours, lb=0, ub=S_FC_max,   name='P_FC')    
    
P_imp     = m.addVars(nHours, lb=0, ub=P_imp_ub,   name="P_imp")               # Imported electricity from the Grid in [W]
P_max_imp = m.addVars(months, lb=0, ub=P_imp_ub,   name="P_max_imp")
P_exp     = m.addVars(nHours, lb=0, ub=P_exp_ub,   name="P_exp")               # Exported electricity to the Grid in [W] # Check upper-bound; #gabriele: replace max(P_ELY) with S_ELY. In general, avoid using max() when defining constraints or design variables.

# Waste Heat recovery varibles
m_cw_ELY = m.addVars(nHours, lb=0, ub=m_cw_ELY_max, name="m_cw")               # Cooling water ELY mass flow
m_cw_FC  = m.addVars(nHours, lb=0, ub=m_cw_FC_max, name="m_cw")   

m_cw_HT  = m.addVars(nHours, lb=0, ub=m_cw_max, name="m_cw_HT")                # Hight Temp water mass flow ELY + FC
m_cw_MT  = m.addVars(nHours, lb=0, ub=m_cw_max, name="m_cw_LT")                # Medium Temp. water mass flow ELY + FC
P_th_HT  = m.addVars(nHours, lb=0, ub=P_th_max, name="P_th_HT")  
P_th_MT  = m.addVars(nHours, lb=0, ub=P_th_max, name="P_th_LT")

S_HP     = m.addVar(lb=0, ub=P_th_max, name="S_HP")
S_HEX    = m.addVar(lb=0, ub=S_HEX_max, name="S_HEX")
# S_HEX  = m.addVar(lb=0, ub=P_th_max, name="S_HEX")  # Uncomment if S_HEX is needed

if include_battery:
    P_ch    = m.addVars(nHours, name="P_ch")                                                    # Power charged to the battery
    P_disch = m.addVars(nHours, name="P_disch")                                                 # Power discharged from the battery
    E_b     = m.addVars(nHours + 1, lb=0, ub=battery_params['C_b_max'], name="E_b")             # Energy in the battery at each time step
    C_b     = m.addVar(lb=battery_params['C_b_min'], ub=battery_params['C_b_max'], name="C_b")  # Battery Capacity

#------------------------------------------------------------------------------
# Additional calculations
#------------------------------------------------------------------------------
# PV calculation
# Sources: De Soto W et al. (2006); Sun V et al. (2020)
def pv_efficiency(irradiance, T_amb):
    eta_PV_ref  = 0.15  # Consider changing back to 0.21 
    T_PV_ref    = 25    # Reference temperature of the pv cell in [°C]
    beta_PV_ref = 0.004 # in [K]
    gamma_PV    = 0.12
    NOCT        = 45
    
    # Calculate T_cell based on ambient Temperature T_amb
    T_cell = T_amb + (NOCT - 20) * irradiance / 800

    # Initialize eff_cell as a zeros array of the same length as irradiance
    eta_cell = np.zeros(len(irradiance))

    # Loop through each irradiance value
    for i in range(len(irradiance)): 
        if irradiance[i] == 0:
            eta_cell[i] = 0
        else:
            eta_cell[i] = eta_PV_ref * (1 - beta_PV_ref * (T_cell[i] - T_PV_ref) + gamma_PV * np.log10(irradiance[i]))
    
    return eta_cell

# Calling the function which calculates the relative efficiency
eta_cell = pv_efficiency(irradiance, T_amb)
# Generate dataframe for plot
df_pv = pd.DataFrame({'irradiance': irradiance,'T_amb': T_amb,'eta_cell': eta_cell})

# # Plot
# plt.figure()
# pc = plt.scatter(df_pv['irradiance'], df_pv['eta_cell'], c=df_pv['T_amb'], cmap='jet')
# plt.colorbar(label='Temperature [C]', ax=plt.gca()) # Adding a colorbar with a label
# pc.set_alpha(0.25) # Setting the alpha for the scatter plot
# plt.grid(alpha=0.5) # Adding grid to the plot with some transparency
# plt.ylim(bottom=0.15)  # # Setting the lower limit for y-axis, adjusted for clarity
# plt.xlabel('Irradiance [W/m²]')
# plt.ylabel('Relative efficiency [-]')
# plt.show()

P_PV = [irradiance[t] * eta_cell[t] * Area_PV  for t in range(nHours)]          # Hourly PV power generation [W]
# P_PV = [irradiance[t] * eta['PV'] * Area_PV  for t in range(nHours)]          # Hourly PV power generation [W]
    
P_PV_peak = 1000 * max(eta_cell) * Area_PV  # Peak power for investment cost [W]
S_PV = P_PV_peak
# Gabriele: do you have a reference for the definition of P_PV_peak you adopted? I normally consider the peak power from Standard Test Conditions (STC), for 
# which irradiance is fixed at 1000 W/m2. I think it is fairer, otherwise you have your peak power depending on sun availability, and given that your unit 
# price is CHF/kW, it would not be fair to have an investiment cost depending on sun availability (how much units you install instead depends on the sun availability)

#------------------------------------------------------------------------------
# Hydrogen mass flow
# Nominal mass flow hydrogen [kg/h]                                          
mdot_H2_nom = (S_ELY * eta["ELY"] * deltat) / HHV                
# Mass flow produced hydrogen [kg/h]                    
mdot_H2 = [(P_ELY[t] * eta['ELY'] * deltat) / HHV for t in range(nHours)] 

# Size (S_C) and operating power (P_C) of the compressor in [W]
S_C = mdot_H2_nom * L_is_C / (eta["C"] * deltat) 
P_C = [mdot_H2[t] * L_is_C / (eta['C']* deltat)  for t in range(nHours)]

#------------------------------------------------------------------------------
# Constraints
#------------------------------------------------------------------------------

for t in range(1, nHours):
    # constraint for the energy balance in time over the storage component
    m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - ((P_FC[t] * deltat) / eta['FC']), name='HESS Balance')  
    
    # constraint for the ELY power not to exceed the storage capacity
    m.addConstr(P_ELY[t]  <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY')       
    
    # constraint for the FC to only use hydrogen available in the storage
    m.addConstr(P_FC[t]   <= E_TANK[t-1] / deltat, name='FC')                   

for t in range(nHours):
    m.addConstr((P_ELY[t]  <= S_ELY),  name="upper_Size_Constraint_ELY")
    m.addConstr((E_TANK[t] <= S_TANK), name= "upper_Size_Constraint_TANK")
    m.addConstr((P_FC[t]   <= S_FC ),  name= "upper_Size_Constraint_FC")
    

m.addConstr(gp.quicksum(P_imp) <= 1 * gp.quicksum(P_demand), name= "GridUse")

# Initializing FC power
m.addConstr(P_FC[0] <= E_TANK[0] / deltat, name= "InitialFC")                                    

# constraint for H2 storage equal at final and last time step (periodicity)
m.addConstr(E_TANK[0] == E_TANK[nHours-1], name='Periodicity') # 08.02: added name to cosntraint


# Overall energy balance: left => consumers | right => generators
if include_battery:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t] + P_ch[t]/battery_params['eff_ch'] + P_exp[t] + P_demand[t] <= P_PV[t] + P_disch[t]*battery_params['eff_disch'] + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance') 
else:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t] + P_exp[t] + P_demand[t] <= P_PV[t] + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance')
    

# Iterate over each month and add constraints based on P_imp[t] for timesteps in that month
for month in months:
    timestep_indices = df_input.index[df_input['MO'] == month].tolist()
    for t in timestep_indices:
        m.addConstr(P_max_imp[month] >= P_imp[t], name=f"max_monthly_constraint_{month}_{t}")
#m.addConstr(P_max_imp == gp.max_(P_imp[t] for t in range (nHours)))

if include_battery:
    # Energy balance for the battery storage
    m.addConstrs(E_b[i + 1] == E_b[i] * battery_params['eff_sd'] + P_ch[i] * battery_params['eff_ch'] - P_disch[i] / battery_params['eff_disch'] for i in range(nHours))
    # Discharge power not exceeding available power
    m.addConstrs(P_disch[i] <= E_b[i] for i in range(nHours))
    # SOC constraints
    m.addConstrs(E_b[i] <= battery_params['SOC_max'] * C_b for i in range(nHours))
    m.addConstrs(E_b[i] >= battery_params['SOC_min'] * C_b for i in range(nHours))
    # Periodicity
    m.addConstr(E_b[nHours] == E_b[0])  # Ensuring energy in the battery at the end matches the start


# Waste Heat Recovery Constraints----------------------------------------------

for t in range(nHours):
    # PEM & FC outlet flow assuming T_HEX and T_in are equivalent for ELY & FC
    m.addConstr(m_cw_ELY[t] == ((1 - eta['ELY']) * P_ELY[t]) / (c_p * (T_HEX - T_in)), "PEM_outlet")
    m.addConstr(m_cw_FC[t] == ((1 - eta['FC']) * P_FC[t]) / (c_p * (T_HEX - T_in)), "FC_outlet")
    
    # Cooling flow requirements
    m.addConstr(m_cw_ELY[t] + m_cw_FC[t] >= m_cw_HT[t] + m_cw_MT[t], "massflowBalanceCoolingWater")
    
    # Low- and High-temperature heat output
    m.addConstr(P_th_MT[t] == m_cw_MT[t] * c_p * (T_out - T_in), "LTheat")
    m.addConstr(P_th_HT[t] == m_cw_HT[t] * c_p * (T_out - T_in), "HTheat")
    
    # Heat exchanger size constraint for medium-temperature side
    m.addConstr(P_th_MT[t] / (U_HEX * T_log) <= S_HEX, "HEXsize")
    
    # Heat pump size constraint for high-temperature side
    m.addConstr(P_th_HT[t] <= S_HP, "HPsize")

# Combined heat exchanger and heat pump constraint
# m.addConstr(P_th_LT + P_th_HT <= eff_th * (P_ELY - (1 - eta['ELY']) * P_ELY), "HEX_HP")

#------------------------------------------------------------------------------
# Run cost function
#------------------------------------------------------------------------------
if include_battery:
    S_BAT = C_b
else:
    S_BAT = 0
    
system_sizes = {'PV': S_PV, 'BAT': S_BAT, 'ELY': S_ELY,'C': S_C,
                'TANK': S_TANK,'FC': S_FC,'HEX': S_HEX,'HP': S_HP}

from cost import totalAnnualCost

[cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, 
cost_maint, cost_WHR] = totalAnnualCost(
                        system_sizes, energy_tariff,
                        UP, maintenance, life, 
                        P_imp, P_max_imp, P_exp, P_th_MT, P_th_HT,
                        cost_imp_el, cost_exp_el, cost_export_heatMT, cost_export_heatHT,
                        df_input, nHours, timeline_choice
                        )

# Total annual costs in [€/y]--------------------------------------------------
cost = cost_inst + cost_op + cost_maint

#------------------------------------------------------------------------------
# Define the objective function: minimal costs
m.setObjective(cost, GRB.MINIMIZE) 

#------------------------------------------------------------------------------
# Solve optimization problem
m.optimize() 

# print (m.display())
# m.printStats()
#------------------------------------------------------------------------------
# Check the optimization status
if m.status == GRB.OPTIMAL:
    print("Optimal solution found.")
elif m.status == GRB.INFEASIBLE:
    print("The model is infeasible.")
elif m.status == GRB.UNBOUNDED:
    print("The model is unbounded.")
elif m.status == GRB.TIME_LIMIT:
    print("Time limit reached.")

# Print the status message for more details => check here: https://www.gurobi.com/documentation/current/refman/optimization_status_codes.html 
print(f"Optimization status: {m.status}")

#------------------------------------------------------------------------------
# Retrieve values of variables for further analysis
#------------------------------------------------------------------------------

# For Gurobi tupledict object
P_ELY  = [P_ELY[t].X  for t in range(nHours)]
P_FC   = [P_FC[t].X   for t in range(nHours)]
E_TANK = [E_TANK[t].X for t in range(nHours)]
P_imp  = [P_imp[t].X  for t in range(nHours)]
P_exp  = [P_exp[t].X  for t in range(nHours)]
P_max_imp = [P_max_imp[month].X for month in months]

if include_battery:
    E_b     = [E_b[t].X for t in range(nHours)]
    P_ch    = [P_ch[t].X for t in range(nHours)]
    P_disch = [P_disch[t].X for t in range(nHours)]
    C_b     = C_b.X
    S_BAT   = S_BAT.X
    

m_cw_ELY = [m_cw_ELY[t].X  for t in range(nHours)]
m_cw_FC  = [m_cw_FC[t].X   for t in range(nHours)]
m_cw_HT  = [m_cw_HT[t].X   for t in range(nHours)]
m_cw_MT  = [m_cw_MT[t].X   for t in range(nHours)]
P_th_HT  = [P_th_HT[t].X   for t in range(nHours)]
P_th_MT  = [P_th_MT[t].X   for t in range(nHours)]

# For lists with Gurobi LinExpr 
P_PV    = [P_PV[t].getValue()    for t in range(nHours)]
P_C     = [P_C[t].getValue()     for t in range(nHours)]
mdot_H2 = [mdot_H2[t].getValue() for t in range(nHours)]

# For Gurobi var object
Area_PV   = Area_PV.X
S_ELY     = S_ELY.X
S_TANK    = S_TANK.X
S_FC      = S_FC.X
S_HEX     = S_HEX.X
S_HP      = S_HP.X

# For Gurobi LinExpr
cost            = cost.getValue()
cost_inst       = cost_inst.getValue()
cost_op         = cost_op.getValue()
cost_maint      = cost_maint.getValue()
cost_elec_imp   = cost_elec_imp.getValue()
cost_elec_exp   = cost_elec_exp.getValue()
cost_grid_usage = cost_grid_usage.getValue()
cost_elec       = cost_elec.getValue()
cost_WHR        = cost_WHR.getValue()

mdot_H2_nom = mdot_H2_nom.getValue()
P_PV_peak   = P_PV_peak.getValue()
S_C         = S_C.getValue()
S_PV        = S_PV.getValue()

all_costs = {
    "installation cost":  cost_inst,
    "imported electricity": cost_elec_imp,
    "exported electricity": cost_elec_exp,
    "grid use fees": cost_grid_usage,
    "total electricity cost": cost_elec,
    "operational cost": cost_op,
    "maintenance cost": cost_maint,
    "revenues WHR": cost_WHR,
    "TAC": cost
    }
#------------------------------------------------------------------------------
S_PV_max = 1000*eta['PV']*Area_PV_max
S_C_max = max(mdot_H2) * L_is_C / (eta["C"] * deltat)
#------------------------------------------------------------------------------
# Calculate the Levelized Cost Of Hydrogen

m_H2_year = sum(mdot_H2) * deltat

# Check if there's any hydrogen generation
if sum(mdot_H2) == 0:
    print("No H2 generation")
else: 
    print("H2 generation")
        

LCOE = (cost) / (sum(P_imp + P_PV)/10**6)  # Calculate LCOE # Look into the literature how
print("LCOE = {:.3f} € / MWh".format(LCOE))

#------------------------------------------------------------------------------
# Display results
print("Area_PV = {:.2f} square meters".format(Area_PV))
#------------------------------------------------------------------------------
# Plot main results------------------------------------------------------------

# Import the plotting module
# from plotting_module_plotly import plot_power_generation, plot_component_sizes, plot_HESS_results, plot_costs_and_prices
from plotting_module import plot_power_generation, plot_component_sizes, plot_HESS_results, plot_battery_operation, plot_costs_and_prices

# Call the plotting functions as needed
# plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours)
plot_component_sizes(S_PV, S_PV_max, S_ELY, S_ELY_max, S_C, S_C_max, S_FC, S_FC_max, S_TANK, S_TANK_max)
plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input)
plot_costs_and_prices(all_costs, df_input)

if include_battery:
    plot_battery_operation(P_demand, P_imp, P_ch, P_disch, E_b, battery_params, C_b, nHours)

fig_power_generation = plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours)
#------------------------------------------------------------------------------
# Export results to excel -----------------------------------------------------

from results_export import export_optimization_results

if include_battery:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'cost_imp_el',
        'P_exp', 'cost_exp_el', 'P_ELY', 'mdot_H2', 'P_C', 'E_TANK',
        'P_FC', 'E_b', 'P_ch', 'P_disch', 'm_cw_ELY', 'm_cw_FC', 'm_cw_HT',
        'm_cw_MT', 'P_th_HT', 'P_th_MT'
        ]
else:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'cost_imp_el',
        'P_exp', 'cost_exp_el', 'P_ELY', 'mdot_H2', 'P_C', 'E_TANK',
        'P_FC', 'm_cw_ELY', 'm_cw_FC', 'm_cw_HT',
        'm_cw_MT', 'P_th_HT', 'P_th_MT'
        ]

results = {name: [] for name in variable_names}

for t in range(nHours):
    results['irradiance'].append(irradiance[t])
    results['P_demand'].append(P_demand[t])
    results['P_PV'].append(P_PV[t])
    results['P_imp'].append(P_imp[t])
    results['cost_imp_el'].append(cost_imp_el[t])
    results['P_exp'].append(P_exp[t])
    results['cost_exp_el'].append(cost_exp_el[t])
    results['P_ELY'].append(P_ELY[t])
    results['mdot_H2'].append(mdot_H2[t])
    results['P_C'].append(P_C[t])
    results['E_TANK'].append(E_TANK[t])
    results['P_FC'].append(P_FC[t])
    if include_battery:
        results['E_b'].append(E_b[t])
        results['P_ch'].append(P_ch[t])
        results['P_disch'].append(P_disch[t])
    results['m_cw_ELY'].append(m_cw_ELY[t])
    results['m_cw_FC'].append(m_cw_FC[t])
    results['m_cw_HT'].append(m_cw_HT[t])
    results['m_cw_MT'].append(m_cw_MT[t])
    results['P_th_HT'].append(P_th_HT[t])
    results['P_th_MT'].append(P_th_MT[t])


# Add optimization status to results
results['status'] = m.status

# Define the path to the results directory
results_directory = export_path

# Export results
export_optimization_results(variable_names, results, results_directory)