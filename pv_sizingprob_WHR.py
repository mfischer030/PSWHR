# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:26:35 2023

@author: fism
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import Model, GRB, LinExpr
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

# path to the functions directory => change in config.py if needed
sys.path.append(function_path) 
#------------------------------------------------------------------------------
# Setting up the model
#------------------------------------------------------------------------------

# Choose energy tariff
energy_tariff = "Green" # Choose from "Green","Blue","Grey"

# Flag to include/exclude battery in the optimization model
include_battery = False

# Choose the "grid-connectivity" of the model: Grid Connectivity Factor
GCF = 20  # [%] of power which can be imported from the grid from the total energy demand

# Capital recovery factor CRF
discountRate = 0.07     # [-] Annual Discount Rate; 2023_Giovanniello | 2021_Marocco

# Define number of breakpoints on PWA approx. of Electrolyser efficiency curve
N_bp = 8

# Import functions and get the efficiency and cost curve coefficients
from efficiencies import pwa_eta_ELY
x_bp_val, y_bp_val, mm_elec, qq_elec = pwa_eta_ELY(N_bp)

from cost_curves import aa_c, bb_c, xx_c, aa_e, bb_e, xx_e

#------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------

# Import data and generate variables using the get_data function---------------
from import_data import get_data
irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand, timeline_choice, season_choice = get_data(input_path, demand_path)


# Defining heat demand variables and convert Wh to kWh by dividing by 1000
heat_zone1_35degC_demand_kWh = df_heat_demand['Heating_Zone1_35degC_W'].values / 1000
heat_zone1_60degC_demand_kWh = df_heat_demand['Heating_Zone1_60degC_W'].values / 1000
heat_zone2_35degC_demand_kWh = df_heat_demand['Heating_Zone2_35degC_W'].values / 1000
heat_zone2_60degC_demand_kWh = df_heat_demand['Heating_Zone2_60degC_W'].values / 1000
heat_zone3_35degC_demand_kWh = df_heat_demand['Heating_Zone3_35degC_W'].values / 1000
heat_zone3_60degC_demand_kWh = df_heat_demand['Heating_Zone3_60degC_W'].values / 1000

# Create a figure and axis for the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each series of data with a label
ax.plot(heat_zone1_35degC_demand_kWh, label='Zone 1, 35°C')
ax.plot(heat_zone1_60degC_demand_kWh, label='Zone 1, 60°C')
ax.plot(heat_zone2_35degC_demand_kWh, label='Zone 2, 35°C')
ax.plot(heat_zone2_60degC_demand_kWh, label='Zone 2, 60°C')
ax.plot(heat_zone3_35degC_demand_kWh, label='Zone 3, 35°C')
ax.plot(heat_zone3_60degC_demand_kWh, label='Zone 3, 60°C')

# Add some plot decorations
ax.set_xlabel('Time')  # Assuming the index represents time
ax.set_ylabel('Demand (kWh)')
ax.set_title('Heating Demand by Zone and Temperature')
ax.legend()

# Show the plot
plt.show()

# # Heat demand data (from NEST building @ empa, year 2022)----------------------
# df_demand_heat = pd.read_excel(heat_path)
# Old plot used for NEST heat demand
# from plotting_module import heat_demand_plot
# heat_demand_plot(heat_35degC_demand,heat_65degC_demand)

#------------------------------------------------------------------------------
# Input parameters
#------------------------------------------------------------------------------
# General parameters
k    = 1.4                               # Ratio cp/cv [-]
R_H2 = 4.1242 * 1000                     # Individual Gas constant H2 [J/kg*K]
M_H2= 2.01568                            # in [g/mol]
HHV  = 39.39 * 3600 * 1000               # Higher heating value of H2 in [J/kg] = 39.39 kWh/kg

# Electricity prices in [EUR/MWh]
cost_imp_el = df_input['price_Eur_MWh'].values # hourly cost to import 
cost_exp_el = df_input['Price_DayAhed'].values # hourly cost to export 

# Revenues from WHR in [€/kWh]: Information on DHN Heat prices can be found here: https://www.preisueberwacher.admin.ch/pue/de/home/themen/infrastruktur/fernwaerme.html
cost_export_heatMT = 0.09   # 0.09
cost_export_heatHT = 0.1189 # 0.1189

# Time parameters
nHours = len(df_input)                        # number of hours simulated
Time   = np.arange(1, nHours + 1)             # time vector
days   = nHours / 24                          # number of days
months = df_input['MO'].unique()              # Unique months in the DataFrame
weeks  = days / 7                             # number of weeks
deltat = 3600                                 # time step (s)

kWh2J = 3600*1000
J2kWh = 1 / (3600*1000)

# Convert nHours to datetime format
if timeline_choice == 'year':
    df_input['ts'] = pd.Timestamp('2019-01-01') + pd.to_timedelta(Time - 1, unit='H')
else:
    if season_choice == 'winter': 
        df_input['ts'] = pd.Timestamp('2019-01-01') + pd.to_timedelta(Time - 1, unit='H')
    elif season_choice == 'summer':
        df_input['ts'] = pd.Timestamp('2019-08-01') + pd.to_timedelta(Time - 1, unit='H')

# Convert df_input['ts'] to a list named 'timeline'
timeline = df_input['ts'].tolist()

# these two coefficients are needed for the big-M constraints
# Tipically, the M value is selected as a large value, but too large values must be avoided to avoid numerical instabilities
# the selection of the M values depend on the problem you are tackling
M_hours     = 2*max(P_demand)    # used for the big-M constraint on the grid_usage
M_threshold = nHours             # used for the big-M for the threshold

# Threshold for low/high grid usage in [h]
threshold_hours = {'year': 3500,'month': 3500/12,'week': 3500/52}
# Small value to avoid numerical issues
epsilon       = 1e-4  

#------------------------------------------------------------------------------
# Efficiencies of components in [-]
#------------------------------------------------------------------------------
"""
PV: from Roxanne
ELY: 2023_van_der_roest (74-79%)
C: from Roxanne | 2020_Pan et al: 0.9
TANK: 2023_Wang et al
FC: 2021_cigolotti Comprehensive Review on FC Technology for Stationary Applications: Electric Efficiency PEMFC: [38,38,37,40] in [%]
"""
if N_bp == 2:
    eta_ELY_nominal = mm_elec + qq_elec # Electrolyser efficiency at nominal power (approx) => from pem_efficiencies.py
else:
    eta_ELY_nominal = mm_elec[-1] + qq_elec[-1]
    
# eta = {'PV': 0.21,'ELY': eta_ELY_nominal,'C': 0.7526,'TANK': 0.95,'FC': 0.5,}              # Realistic
#eta = {'PV': 0.21,'ELY': 0.8,'C': 0.8763,'TANK': 0.975,'FC': 0.75}
eta = {'PV': 0.21,'ELY': 1,'C': 1,'TANK': 1,'FC': 1}                                         # Optimal      

#------------------------------------------------------------------------------
# Unit prices of components / capital costs in [€/W] 
#------------------------------------------------------------------------------
"""
PV:   in [€/W] | 2023_Tay Son Le: 881 USD/year | 2018_Gabrielli 300€/m2 
BAT:  in [€/kWh]
ELY:  in [€/W] | 2023_IEA-GlobalH2Review: PEM - 2kUSD/kW - reduction to 600 USD/kW
C:    in [€/W] | 2020_Pan et al: 1228 (¥/kW)
TANK: in [EUR/J], 1644 [EUR/kgH2] | 2018_Gabrielli: {20.7;13.6;10.9} in [€/kWh] | 2023_A Review on the Cost Analysis of H2 Gas Storage Tanks for FCV: 2020:9.34 €/kWh 2025: 8.40€/kWh Ultimate: 7.47€/kWh
FC:   in [€/W] (Wang et. al 2023 - 2kUSD/kW) | 2018_Gabrielli {2160;1680;1320} [€/kW]
HEX:  in [€/m2] from Roxanne: 77.79 €/m2 + Fixed_HEX = 5291.9 (Fixed cost for heat exchanger [EUR]) | From Christian: 1782CHF/5.2m2 = 342.69 CHF/m2
"""


# UP = {'PV': 0.8,'BAT': 0.5/3600,'ELY': 1.7,'C': 0.0076, 'TANK': 9.34/(3.6*10**6),'FC': 2,  'HP': 0.576, 'HEX': 342.69}      # Realistic
# UP = {'PV': 0.8,'BAT': 0.5/3600,'ELY': 1.128,'C': 0.00506,'TANK': 8.4/(3.6*10**6), 'FC': 1.3,'HP': 0.238, 'HEX': 221}     # Mid-Term
UP  =  {'PV': 0.8,'BAT': 0.5/3600,'ELY': 0.556,'C': 0.0038, 'TANK': 7.47/(3.6*10**6),'FC': 0.6,'HP': 0.200, 'HEX': 100}   # Optimal


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

# Daily energy demand in [Wh/day] - P_demand is in [W] integrate over 8760 hours => [Wh]
E_demand_day = sum(P_demand)/days # Gabriele: this is [Wh], correct, Maxime: yes

# Maximal daily power demand in [W]
P_peak_max = np.max(P_demand)

# Battery Parameters
bat_params = {
    'eff_ch': 0.95,         # Battery charging efficiency in [-]
    'eff_disch': 0.95,      # Battery discharging efficiency in [-] 
    'eff_sd': 0.99,         # Battery self-discharge efficiency in [-]
    'C_b_max': 1000*kWh2J,  # Max battery capacity in [J]
    'C_b_min': 0,           # Min battery capacity in [J] - typically set to 0
    'SOC_max': 0.8,         # Max state of charge of the battery to increase lifetime
    'SOC_min': 0.2          # Min state of charge of the battery to increase lifetime
}

#------------------------------------------------------------------------------
# Defining maximal SIZES of the components
#------------------------------------------------------------------------------

# Area_PV_max  = P_peak_max / (eta['PV'] * np.mean(irradiance)) # Maximum PV area [m2]
Area_PV_max  = 5000 # Maximum PV area [m2] => how to set, quantify the maximum pv area?

# Electrolyser max and min nominal power (W)   
S_ELY_max = 500*1000   # Maximal Electrolyzer size in [W]
S_ELY_min = 1           # Min. size ELY where problem is feasible [W] - from Rox

# Calculating the spezific work of the compresso, from Minutillo et al. 2021
# (Analyzing the levelized cost of hydrogen eq 1+2) => from Roxanne 
T_in_H2 = 65 + 273.15                       # H2 temperature (=T_cat=T_an) [K]  
p_out   = 30                                # Compressor outlet pressure [bar] = H2 storage pressure 
p_in    = 30                                # Compressor inlet pressure [bar]  PEM electrolyzer at Empa works at 30 bar 
L_is_C  = (k/(k-1)) * R_H2 * T_in_H2 * (((p_out/p_in)**((k-1)/k)) - 1) # Specific work compressor [J/kg] 

# Maximal TANK energy capacity (J) => 4 days storage capacity
S_TANK_max = E_demand_day * 30 * 3600                     # E_demand_day in [Wh]  
S_TANK_H2_max = S_TANK_max / HHV                          # equivalent in kg_H2

# Maximal FC size as the maximal power demand divided by eff in [W]
#S_FC_max = P_peak_max / eta["FC"]
S_FC_max = 500*1000  
"""
Range PEMFC = [10W;1MW] from 2021_cigolotti Comprehensive Review on Fuel Cell 
Technology for Stationary Applications
"""
# Waste Heat Recovery ---------------------------------------------------------

T_in  = 57    # Inlet temperature cooling water to HEX in [°C]
T_out = 62    # Temperature cooling water to applications [°C]
T_HEX = 64    # Outlet temperature cooling water from HEX [°C]
c_p   = 4186  # Specific heat capacity cooling water [J/kgK]

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
U_HEX     = 2000      # Overall heat transfer coeff HEX [W/m2*K], swedish thesis

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
#------------------------------------------------------------------------------

m = Model()

#------------------------------------------------------------------------------
# Design variables
#------------------------------------------------------------------------------

# Sizing Variables
Area_PV = m.addVar(lb=0,         ub=Area_PV_max, name='Area_PV')
S_ELY   = m.addVar(lb=S_ELY_min, ub=S_ELY_max,   name='S_ELY')
S_TANK  = m.addVar(lb=0,         ub=S_TANK_max,  name='S_TANK')
S_FC    = m.addVar(lb=0,         ub=S_FC_max,    name='S_FC')

S_HP     = m.addVar(lb=0, ub=P_th_max, name="S_HP")
S_HEX    = m.addVar(lb=0, ub=S_HEX_max, name="S_HEX")
# S_HEX  = m.addVar(lb=0, ub=P_th_max, name="S_HEX")  # Uncomment if S_HEX is needed

# Operation Variables 
P_ELY     = m.addVars(nHours, lb=0, ub=S_ELY_max, name='P_ELY')
P_ELY_PWA = m.addVars(nHours, lb=0, ub=S_ELY_max, name='P_ELY')                  # Output from the Electrolyser => P_ELY[t] * eta[ELY][t]
P_ElyOn   = m.addVars(nHours, lb=0, ub=S_ELY_max, name="P_ElyOn")
ElyOn     = m.addVars(nHours, lb=0, ub=1,vtype=GRB.INTEGER, name="ElyOn")        # Binary Variable for ON-OFF condition of the Electrolyser

E_TANK  = m.addVars(nHours, lb=0, ub=S_TANK_max, name='E_TANK')

# Parameters provided
P_bp     = 47.97 * 1000         # Breakpoint power in Watts (converted from kW)
i_bp     = 123.05               # Breakpoint current in A
P_min_fc = 30 * 1000            # Minimum power in Watts (converted from kW)
P_max_fc = 1000 * 1000          # Maximum power in Watts (converted from kW)

P_FC = m.addVars(nHours, lb=0, ub=S_FC_max, name='P_FC')         # Fuel Cell Output Power in [W]

# P_FC_PWA = m.addVars(nHours, lb=0, ub=S_FC_max, name='P_FC_PWA')         # Fuel Cell Output Power in [W]
# u_FC     = m.addVars(nHours, lb=0, ub=S_FC_max, name='P_FC')             # Fuel Cell Input Power in [W]
# i_FC     = m.addVars(nHours, name='i_fc')
# Vdot_FC  = m.addVars(nHours, name='v_fc')                                # H2 Volume Flow in [Nm3/s]  
# z1       = m.addVars(nHours, vtype=GRB.BINARY, name='z1')                # Binary Variable for PWA (slope 1)
# z2       = m.addVars(nHours, vtype=GRB.BINARY, name='z2')                # Binary Variable for PWA (slope 2)
    
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


if include_battery:
    P_ch    = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="P_ch")              # Power charged to the battery
    P_ds    = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="P_ds")              # Power discharged from the battery
    E_b     = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="E_b")               # Energy in the battery at each time step
    C_b     = m.addVar(lb=bat_params['C_b_min'], ub=bat_params['C_b_max'], name="C_b")   # Battery Capacity
    

grid_usage = m.addVars(nHours, vtype=GRB.BINARY, name="grid_usage")  # Binary variable for grid usage indicator, 1 if used, i.e. P_imp>0, 0 otherwise
h_usage    = m.addVar(lb=0, ub=nHours, name="h_usage")               # Total hours of grid usage
high_usage = m.addVar(vtype=GRB.BINARY, name="high_usage")           # Binary variable indicating whether grid usage is above the threshold, 1 if above, 0 if below


#------------------------------------------------------------------------------
# Additional calculations
#------------------------------------------------------------------------------

# PV calculation---------------------------------------------------------------

# Load the pv_efficiency function from efficiencies.py
from efficiencies import pv_efficiency

# Calling the function which calculates the relative efficiency
eta_cell = pv_efficiency(irradiance, T_amb)

# Generate dataframe for plot
df_pv = pd.DataFrame({'irradiance': irradiance,'T_amb': T_amb,'eta_cell': eta_cell})

P_PV = [irradiance[t] * eta_cell[t] * Area_PV  for t in range(nHours)]          # Hourly PV power generation [W]
# P_PV = [irradiance[t] * eta['PV'] * Area_PV  for t in range(nHours)]          # Hourly PV power generation [W]
    
P_PV_peak = 1000 * max(eta_cell) * Area_PV  # Peak power for investment cost [W]
S_PV = P_PV_peak

# Gabriele: do you have a reference for the definition of P_PV_peak you adopted? I normally consider the peak power from Standard Test Conditions (STC), for 
# which irradiance is fixed at 1000 W/m2. I think it is fairer, otherwise you have your peak power depending on sun availability, and given that your unit 
# price is CHF/kW, it would not be fair to have an investiment cost depending on sun availability (how much units you install instead depends on the sun availability)

# Hydrogen mass flow-----------------------------------------------------------

# Nominal mass flow hydrogen [kg/h]                                          
mdot_H2_nom = S_ELY     * eta["ELY"] / HHV 
mdot_H2_max = S_ELY_max * eta['ELY'] / HHV               

# Mass flow produced hydrogen [kg/h]                    
# mdot_H2 = [P_ELY[t] * eta['ELY']  / HHV for t in range(nHours)] 
mdot_H2   = [P_ELY_PWA[t] / HHV for t in range(nHours)] 

# Size (S_C) and operating power (P_C) of the compressor in [W]
S_C     = mdot_H2_nom * L_is_C / eta["C"]  
S_C_max = mdot_H2_max * L_is_C / eta["C"] 
P_C     = [mdot_H2[t] * L_is_C / eta['C'] for t in range(nHours)]

#------------------------------------------------------------------------------
# Constraints
#------------------------------------------------------------------------------

for t in range(1, nHours):
    # Energy Balance for HESS
    # m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - ((P_FC[t] * deltat) / eta['FC']), name='HESS Balance')
    m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY_PWA[t] * deltat - ((P_FC[t] * deltat) / eta['FC']), name='HESS Balance')
    # m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY_PWA[t] * deltat - P_FC_PWA[t] * deltat, name='HESS Balance')  
    
    # constraint for the ELY power not to exceed the storage capacity
    # m.addConstr(P_ELY[t]  <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY') 
    m.addConstr(P_ELY_PWA[t]  <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY')       
    
    # constraint for the FC to only use hydrogen available in the storage
    m.addConstr(P_FC[t]   <= E_TANK[t-1] / deltat, name='FC')                   

for t in range(nHours):
    m.addConstr((P_ELY[t]  <= S_ELY),  name="upper_Size_Constraint_ELY")
    m.addConstr((E_TANK[t] <= S_TANK), name= "upper_Size_Constraint_TANK")
    m.addConstr((P_FC[t]   <= S_FC ),  name= "upper_Size_Constraint_FC")

# Initializing FC power
m.addConstr(P_FC[0] <= E_TANK[0] / deltat, name= "InitialFC")                                    

# constraint for H2 storage equal at final and last time step (periodicity)
m.addConstr(E_TANK[0] == E_TANK[nHours-1], name='Periodicity_HESS') # 08.02: added name to cosntraint

# Overall energy balance: left => consumers | right => generators
if include_battery:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t] + P_ch[t] + P_exp[t] + P_demand[t] <= P_PV[t] + P_ds[t] + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance') 
else:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t]           + P_exp[t] + P_demand[t] <= P_PV[t]           + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance')

# Constraint for the Grid Connectivity 
for t in range(nHours):
    # m.addConstr(P_imp[t] <= GCF/100 * P_demand[t], name= "GridUse")
    m.addConstr(gp.quicksum(P_imp) <= GCF/100 * gp.quicksum(P_demand), name= "GridUse")

# Monthly import peaks for cost function
for month in months:
    timestep_indices = df_input.index[df_input['MO'] == month].tolist()
    for t in timestep_indices:
        m.addConstr(P_max_imp[month] >= P_imp[t], name=f"max_monthly_constraint_{month}_{t}")

# Battery Constraints ---------------------------------------------------------
if include_battery:
    # Energy balance for the battery storage
    for t in range(1, nHours):
        m.addConstr(E_b[t] == E_b[t-1] * bat_params['eff_sd'] + P_ch[t] * bat_params['eff_ch'] * deltat - ((P_ds[t] * deltat) / bat_params['eff_disch']), name='Battery energy balance')
    
    for t in range(nHours):
        # Discharge power not exceeding available power
        m.addConstr(P_ds[t] <= E_b[t])
    
        # SOC constraints
        m.addConstr(E_b[t] <= bat_params['SOC_max'] * C_b)
        # m.addConstr(E_b[t] >= bat_params['SOC_min'] * C_b)
    
    # Periodicity = Ensuring energy in the battery at the end matches the start
    m.addConstr(E_b[0] == E_b[nHours-1], name='Periodicity_Battery')

# Waste Heat Recovery Constraints----------------------------------------------

for t in range(nHours):
    # PEM & FC outlet flow assuming T_HEX and T_in are equivalent for ELY & FC
    m.addConstr(m_cw_ELY[t] == ((1 - eta['ELY']) * P_ELY[t]) / (c_p * (T_HEX - T_in)), "PEM_outlet")
    m.addConstr(m_cw_FC[t]  == ((1 - eta['FC'])  * P_FC[t])  / (c_p * (T_HEX - T_in)), "FC_outlet")
    
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
# m.addConstr(P_th_T + P_th_HT <= eff_th * (P_ELY - (1 - eta['ELY']) * P_ELY), "HEX_HP")

# BIG-M constraint for Benutzungsdauer = grid usage time-----------------------
# we use the big-M approach to define the binary variable grid_usage
for i in range(nHours):
    m.addConstr(P_imp[i] <= M_hours * grid_usage[i])                 # here we constraint grid_usage to take value 1 if P_imp is larger than zero. If P_imp is positive, this constraint is satisfied only if grid_usage is 1
    m.addConstr(P_imp[i] >= epsilon - M_hours * (1 - grid_usage[i])) # here we constraint grid_usage to 0 if P_imp is zero. 
# Here I choose M_hours in the order of the peak demand

# here we add the constraints needed for the grid usage tarif selection
m.addConstr(h_usage == gp.quicksum(grid_usage[i] for i in range(nHours))) # not sure this has to be a constraint, could be a normal calculation maybe 

# defining the value of the binary variable high-usage based on the hours of operation, h_usage
m.addConstr(h_usage - threshold_hours[timeline_choice] <= M_threshold * high_usage)# if h_usage is larger than threshold, the inequality is respected only if high_usage takes value 1. 
# here I choose the M_threshold in the order of the hours in the year

# PWA approximation of the Electrolyser----------------------------------------

for t in range(nHours):
    # Constraint: P_e should be less than or equal to P_ElyOn
    m.addConstr(P_ELY[t] <= P_ElyOn[t], "MaxPowerEly")

    # Constraint: P_e should be at least 20% of P_ElyOn
    m.addConstr(P_ELY[t] >= 0.2 * P_ElyOn[t], "MinPowerEly")

    # Constraint: P_ElyOn not to exceed the maximum power when the electrolyser is on (S_e_max) 
    m.addConstr(P_ElyOn[t] <= S_ELY_max * ElyOn[t], "PowerEly1")

    # Constraint: P_ElyOn should not exceed S_e
    m.addConstr(P_ElyOn[t] <= S_ELY, "PowerEly3")

    # Constraint: When the electrolyser is off (ElyOn = 0), P_ElyOn should be at least S_e - S_e_max. When ElyOn = 1, this constraint has no effect because S_e - S_e_max * (1 - ElyOn) equals S_e - S_e_max * 0 = S_e
    m.addConstr(P_ElyOn[t] >= S_ELY - S_ELY_max * (1 - ElyOn[t]), f"PowerEly4_{i}")
    
    # Piecewise linear approximation constraints
    if N_bp == 1:
        m.addConstr(P_ELY_PWA[t] == eta_ELY_nominal * P_ELY[t], "Efficiency Electrolyser N_bp=1")
    elif N_bp == 2:
        m.addConstr(P_ELY_PWA[t] <= mm_elec * P_ELY[t] + qq_elec * P_ElyOn[t],  "PWA_Q1")
    elif N_bp==3:
        m.addConstr(P_ELY_PWA[t] <= mm_elec[0] * P_ELY[t] + qq_elec[0] * P_ElyOn[t],  "PWA_Q1")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[1]  * P_ELY[t] + qq_elec[1] * P_ElyOn[t], "PWA_Q2")
    elif N_bp==4:
        m.addConstr(P_ELY_PWA[t] <= mm_elec[0] * P_ELY[t] + qq_elec[0] * P_ElyOn[t],  "PWA_Q1")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[1]  * P_ELY[t] + qq_elec[1] * P_ElyOn[t], "PWA_Q2")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[2]  * P_ELY[t] + qq_elec[2] * P_ElyOn[t], "PWA_Q3")
    elif N_bp==8:
        m.addConstr(P_ELY_PWA[t] <= mm_elec[0]  * P_ELY[t] + qq_elec[0] * P_ElyOn[t],  "PWA_Q1")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[1]  * P_ELY[t] + qq_elec[1] * P_ElyOn[t], "PWA_Q2")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[2]  * P_ELY[t] + qq_elec[2] * P_ElyOn[t], "PWA_Q3")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[3]  * P_ELY[t] + qq_elec[3] * P_ElyOn[t],  "PWA_Q4")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[4]  * P_ELY[t] + qq_elec[4] * P_ElyOn[t],  "PWA_Q5")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[5]  * P_ELY[t] + qq_elec[5] * P_ElyOn[t],  "PWA_Q6")
        m.addConstr(P_ELY_PWA[t] <= mm_elec[6]  * P_ELY[t] + qq_elec[6] * P_ElyOn[t],  "PWA_Q7")
        
    # m.addConstr(P_ELY_PWA[t] <= mm_elec[7]  * P_ELY[t] + qq_elec[7] * P_ElyOn[t],  "PWA_Q8")
    # m.addConstr(P_ELY_PWA[t] <= mm_elec[8]  * P_ELY[t] + qq_elec[8] * P_ElyOn[t],  "PWA_Q9")
    # m.addConstr(P_ELY_PWA[t] <= mm_elec[9]  * P_ELY[t] + qq_elec[9] * P_ElyOn[t],  "PWA_Q10")
    # m.addConstr(P_ELY_PWA[t] <= mm_elec[10] * P_ELY[t] + qq_elec[10] * P_ElyOn[t], "PWA_Q11")

# PWA approximation of the Fuel Cell-------------------------------------------

s1       = 2.56 * 1000          # Slope below breakpoint in [1/V]
s2       = 3.31 * 1000          # Slope above breakpoint in [1/V]
c        = 0.21                 # Volume flow conversion factor in [Nm3/C]
d        = 0.96                 # Efficiency coefficient [-]

# Constraints
# for t in range(nHours):
#     # Enforce the binary variables for the segments
#     m.addConstr(z1[t] + z2[t] == 1, name=f'piecewise_selection_{t}')

#     # First segment constraints
#     m.addConstr(i_FC[t] <= s1 * u_FC[t] + (1 - z1[t]) * S_FC_max, name=f'first_segment_upper_{t}')
#     m.addConstr(i_FC[t] >= s1 * u_FC[t] - (1 - z1[t]) * S_FC_max, name=f'first_segment_lower_{t}')
#     m.addConstr(u_FC[t] <= P_bp + (1 - z1[t]) * S_FC_max, name=f'first_segment_range_{t}')

#     # Second segment constraints
#     m.addConstr(i_FC[t] <= s2 * (u_FC[t] - P_bp) + i_bp + (1 - z2[t]) * S_FC_max, name=f'second_segment_upper_{t}')
#     m.addConstr(i_FC[t] >= s2 * (u_FC[t] - P_bp) + i_bp - (1 - z2[t]) * S_FC_max, name=f'second_segment_lower_{t}')
#     m.addConstr(u_FC[t]  >= P_bp - (1 - z2[t]) * S_FC_max, name=f'second_segment_range_{t}')

#     # Calculate the H2 volume flow 'Vdot_FC' based on the current 'i_FC'
#     m.addConstr(Vdot_FC[t] == c * i_FC[t], name=f'H2volume_flow_{t}')

#     # Calculate the PWA power output 'P_FC_PWA' based on the power input 'u_FC'
#     m.addConstr(P_FC_PWA[t] == d * u_FC[t], name=f'power_output_{t}')

# Implement cost curves: compressor and electrolyser --------------------------


# cost_C_max = aa_c[1]*S_C_max + bb_c[1]
# cost_C     = m.addVar(lb=0, ub=cost_C_max, name="cost_c")
# y_c1       = m.addVar(vtype=GRB.BINARY,    name="y_c1")
# y_c2       = m.addVar(vtype=GRB.BINARY,    name="y_c2")
# S_cy1      = m.addVar(lb=0, ub=S_C_max,    name="S_cy1")
# S_cy2      = m.addVar(lb=0, ub=S_C_max,    name="S_cy2")
    
# # Comp constraints
# m.addConstr(S_cy1 <= S_C_max * y_c1,             "Comp_cy11")
# m.addConstr(S_cy1 <= S_C,                        "Comp_cy13")
# m.addConstr(S_cy1 >= S_C - S_C_max * (1 - y_c1), "Comp_cy14")
    
# m.addConstr(S_cy2 <= S_C_max * y_c2,             "Comp_cy21")
# m.addConstr(S_cy2 <= S_C,                        "Comp_cy23")
# m.addConstr(S_cy2 >= S_C - S_C_max * (1 - y_c2), "Comp_cy24")
    
# m.addConstr(cost_C == aa_c[0] * S_cy1 + bb_c[0] * y_c1 + aa_c[1] * S_cy2 + bb_c[1] * y_c2, "Comp_cost")
# m.addConstr(y_c1 + y_c2 == 1, "Comp_cost2")
# m.addConstr(S_C <= xx_c[0] * y_c1 + xx_c[1] * y_c2, "Comp_cost4")


# Define the breakpoints for the piecewise approximation of the cost curves
x_bp_ELY = [250000, 400000, 550000]
y_bp_ELY = [347388, 495479, 632884]

x_bp_C = [5000, 12500, 20000]
y_bp_C = [112682, 192791, 253934]

# Define the 'S_C' variables which represent the sizes of C
S_C = m.addVar(vtype=GRB.CONTINUOUS, name="S_C")

# Define binary variables for selecting the piecewise segments for ELY and C
y_ELY = m.addVars(len(x_bp_ELY) + 1, vtype=GRB.BINARY, name="y_ELY")
y_C   = m.addVars(len(x_bp_C) + 1, vtype=GRB.BINARY, name="y_C")

# Ensure only one y variable for ELY and C can be 1 at the same time
m.addConstr(sum(y_ELY[i] for i in range(len(y_ELY))) == 1, "UniqueSegment_ELY")
m.addConstr(sum(y_C[i]   for i in range(len(y_C)))   == 1, "UniqueSegment_C")

# Add the piecewise linear cost functions using the breakpoints
cost_ELY = m.addVar(vtype=GRB.CONTINUOUS, name="cost_ELY")
m.addGenConstrPWL(S_ELY, cost_ELY, x_bp_ELY, y_bp_ELY, "PWL_ELY")

# # Create binary variable for the condition S_C > 0
# is_S_C_positive = m.addVar(vtype=GRB.BINARY, name="is_S_C_positive")

# # Set is_S_C_positive to 1 if S_c is greater than a very small number (effectively 0)
# m.addGenConstrIndicator(is_S_C_positive, True, S_C >= 1e-4, name="Indicator_S_C_positive")

# # Add the piecewise linear cost functions using the breakpoints for C
# cost_C = m.addVar(vtype=GRB.CONTINUOUS, name="cost_C")

# # Link the cost_C to the condition that is_S_C_positive is True
# m.addGenConstrIndicator(is_S_C_positive, True, m.addGenConstrPWL(S_C, cost_C, x_bp_C, y_bp_C, "PWL_C"), name="Cost_C_PWL")

# # Add a constraint to set cost_C to 0 when S_C <= 0
# m.addConstr((1 - is_S_C_positive) * cost_C == 0, "Cost_C_Zero")

#------------------------------------------------------------------------------
# Run cost function
#------------------------------------------------------------------------------
if include_battery:
    S_BAT = C_b
else:
    S_BAT = 0
    
system_sizes = {'PV': S_PV,
                'BAT': S_BAT,
                'HESS':{'ELY': S_ELY,'C': S_C,'TANK': S_TANK,'FC': S_FC},
                'WHR':{'HEX': S_HEX,'HP': S_HP}
                }

from cost import totalAnnualCost

# [cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, 
#   cost_maint, cost_WHR, electricity_prices, costs_with_annuity] = totalAnnualCost(
#                         system_sizes, energy_tariff, discountRate,
#                         UP, maintenance, life, 
#                         P_imp, P_max_imp, P_exp, P_th_MT, P_th_HT,
#                         cost_imp_el, cost_exp_el, cost_export_heatMT, cost_export_heatHT,
#                         m, high_usage,
#                         df_input, nHours, timeline_choice,
#                         cost_ELY, cost_C
#                         )
     
[cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, 
  cost_maint, cost_WHR, electricity_prices] = totalAnnualCost(
                        system_sizes, energy_tariff, discountRate,
                        UP, maintenance, life, 
                        P_imp, P_max_imp, P_exp, P_th_MT, P_th_HT,
                        cost_imp_el, cost_exp_el, cost_export_heatMT, cost_export_heatHT,
                        m, high_usage,
                        df_input, nHours, timeline_choice,
                        cost_ELY, #)cost_C
                        )

# Total annual costs in [€/y]--------------------------------------------------

# Defining a factor to adjust overall cost to choosen timeline
if timeline_choice == 'week':
    multiplier = 52                               # Assuming 52 weeks in a year
elif timeline_choice == 'month':
    multiplier = 12                               # 12 months in a year
else:
    multiplier = 1                                # Default to yearly calculation with no multiplication needed

cost = (cost_inst + cost_op + cost_maint)/multiplier

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
P_ELY   = [P_ELY[t].X   for t in range(nHours)]
P_ElyOn = [P_ElyOn[t].X for t in range(nHours)]
P_ELY_PWA = [P_ELY_PWA[t].X for t in range(nHours)]
ElyOn   = [ElyOn[t].X   for t in range(nHours)]

P_FC = [P_FC[t].X for t in range(nHours)]

# P_FC_PWA = [P_FC_PWA[t].X for t in range(nHours)]
# u_FC     = [u_FC[t].X     for t in range(nHours)]
# i_FC     = [i_FC[t].X     for t in range(nHours)]
# Vdot_FC  = [Vdot_FC[t].X  for t in range(nHours)]
# z1       = [z1[t].X       for t in range(nHours)]
# z2       = [z2[t].X       for t in range(nHours)]


E_TANK  = [E_TANK[t].X  for t in range(nHours)]
P_imp   = [P_imp[t].X   for t in range(nHours)]
P_exp   = [P_exp[t].X   for t in range(nHours)]
P_max_imp = [P_max_imp[month].X for month in months]

if include_battery:
    E_b     = [E_b[t].X for t in range(nHours)]
    P_ch    = [P_ch[t].X for t in range(nHours)]
    P_ds = [P_ds[t].X for t in range(nHours)]
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
Area_PV    = Area_PV.X
S_ELY      = S_ELY.X
S_C        = S_C.X
S_TANK     = S_TANK.X
S_FC       = S_FC.X
S_HEX      = S_HEX.X
S_HP       = S_HP.X
high_usage = high_usage.X
h_usage    = h_usage.X

# cost_C   = cost_C.X
cost_ELY = cost_ELY.X

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
# S_C         = S_C.getValue()
S_PV        = S_PV.getValue()



def retrieve_values(costs):
    """Recursively retrieve values from a nested dictionary of Gurobi LinExpr objects and floats."""
    retrieved_costs = {}
    for key, value in costs.items():
        if isinstance(value, dict):
            # Recurse into the sub-dictionary
            retrieved_costs[key] = retrieve_values(value)
        elif hasattr(value, 'getValue'):
            # Check if the object has the 'getValue' method typical of Gurobi LinExpr objects
            retrieved_costs[key] = value.getValue()
        else:
            # Handle floats and other types that do not need the getValue method
            retrieved_costs[key] = value
    return retrieved_costs

# Assuming costs_with_annuity is defined and some entries might be Gurobi LinExpr objects or floats
# costs_with_annuity = retrieve_values(costs_with_annuity)

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
# Post processing
#------------------------------------------------------------------------------
#Calculate the Benutzungsdauer
BD = sum(P_imp)/np.mean(P_max_imp)


S_PV_max = 1000*eta['PV']*Area_PV_max
S_C_max = max(mdot_H2) * L_is_C / (eta["C"] * deltat)

# Convert the battery capacity from Joules to Kilowatthours
if include_battery:
    C_b_kWh = C_b * J2kWh 

# Calculate the electricity price from BKW ------------------------------------
energy_price = electricity_prices['Electricity prices [Rp./kWh]']
grid_fee     = electricity_prices['Grid use cost [Rp./kWh]']
export_price = electricity_prices['Elecricity export price [Rp.kWh]']

# Import and Export prices in CHF/MWh
electricity_price_imp = [(ip + gf) * 10 for ip, gf in zip(energy_price, grid_fee)]
electricity_price_exp = [ep        * 10 for ep in export_price]
#------------------------------------------------------------------------------
# Function to calculate the volume of hydrogen gas based on temperature and pressure
# Important notice: the coefficients A-E from the paper might not be the right ones this study.
def hydrogen_tank_volume(p_out, T_amb, E_TANK, HHV, R_H2, M_H2, nHours):
    """
    Calculate the volume of hydrogen gas based on the ideal gas law, considering the compressibility factor.
    :param p_out: Electrolyser outlet pressure in bar
    :param T_amb: Ambient Temperature in °C (list of temperatures for each hour)
    :param mdot_H2: Mass flow rate of hydrogen for each hour
    :param R_H2: Ideal gas constant for hydrogen
    :param M_H2: Molar mass of hydrogen
    :param nHours: Number of hours (length of the time series data)
    :return: Volume in m3? (list of volumes for each hour)
    """
    # Constants for the Z factor calculation
    A = 4.93482 * 10**(-5)
    B = 2.04036
    C = 8.15334 * 10
    D = -6.5561 * 10**4
    E = 4.56516 * 10**6
    
    V_Tank_H2 = [0] * nHours  # Initialize hydrogen storage volume list
    
    for t in range(nHours):
        T_K = T_amb[t] + 273.15                 # Convert temperature from °C to Kelvin
        n   = (E_TANK[t] / HHV) / (M_H2/1000)   # Calculate the amount of hydrogen in mol
        if p_out < 500: #!!!!!!!!! CHANGE TO 13 bar which is the normal value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Z = 1  # Use Z = 1 for p_out < 13 bar
        else:
            # Calculate Z for p_out >= 13 without creating a list
            Z = 1 + (p_out*100000) * (A + B * T_K**(-1) + C * T_K**(-2) + D * T_K**(-3) + E * T_K**(-4))
        
        # Calculate the volume using the modified ideal gas law: PV = ZnRT
        V_Tank_H2[t] = Z * n * R_H2 * T_K / (p_out * 100000)
    
    return V_Tank_H2, Z

if sum(mdot_H2) > 0:
    V_Tank_H2, Z = hydrogen_tank_volume(p_out, T_amb, E_TANK, HHV, R_H2, M_H2, nHours)

#------------------------------------------------------------------------------
# Calculate the Levelized Cost Of Hydrogen
#------------------------------------------------------------------------------

m_H2_year = sum(mdot_H2) * deltat

# Check if there's any hydrogen generation
if sum(mdot_H2) == 0:
    print("No H2 generation")
else: 
    print("H2 generation")

LCOE = (cost) / (sum(P_imp + P_PV)/10**6)  # Calculate LCOE # Look into the literature how
VALCOE = cost / (sum(P_demand)/10**6)
#------------------------------------------------------------------------------
# Display results
#------------------------------------------------------------------------------
simulation_report = {
    "Choosen Timeline": timeline_choice,
    "Chosen energy tariff": energy_tariff,
    "Include battery in the simulation": include_battery,
    "Grid Connectivity Factor (GCF) [%]": GCF,
    "Components efficiencies": eta,
    "Components unit prices": UP,
    "Maximal Area covered by PV panels [m2]": Area_PV_max,
    "H2 Storage Pressure [bar]": p_out
}

# Convert dictionary to DataFrame for better visualization
df_simulation_report = pd.DataFrame(list(simulation_report.items()), columns=['Parameter', 'Value'])
# print(df_simulation_report)

# Create a figure and an axes to display the table with adjusted settings
fig, ax = plt.subplots(figsize=(15, 4))  # Increased figure size for better content fit
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=df_simulation_report.values, colLabels=df_simulation_report.columns, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(16)  # Slightly reduced font size
table.scale(1, 2)  # Adjusted cell scaling for better visibility and to prevent overlap
plt.title('Simulation Report')
plt.show()



# Generate a dictionnary with the results
if include_battery:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
        'P_exp', 'electricity_price_exp', 'P_ELY', 'P_ElyOn', 'ElyOn', 
        'P_ELY_PWA', 'mdot_H2', 'P_C', 'E_TANK','P_FC', 'E_b', 'P_ch', 'P_ds', 
        'm_cw_ELY', 'm_cw_FC', 'm_cw_HT','m_cw_MT', 'P_th_HT', 'P_th_MT', 
        'heat_zone3_35degC_demand_kWh','heat_zone3_60degC_demand_kWh'
        ]
else:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
        'P_exp', 'electricity_price_exp', 'P_ELY', 'P_ElyOn', 'ElyOn', 
        'P_ELY_PWA', 'mdot_H2', 'P_C', 'E_TANK',
        'P_FC', 'm_cw_ELY', 'm_cw_FC', 'm_cw_HT',
        'm_cw_MT', 'P_th_HT', 'P_th_MT', 'heat_zone3_35degC_demand_kWh',
        'heat_zone3_60degC_demand_kWh'
        ]

results = {name: [] for name in variable_names}

results['ts'] = timeline

for t in range(nHours):
    results['irradiance'].append(irradiance[t])
    results['P_demand'].append(P_demand[t])
    results['P_PV'].append(P_PV[t])
    results['P_imp'].append(P_imp[t])
    results['electricity_price_imp'].append(electricity_price_imp[t])
    results['P_exp'].append(P_exp[t])
    results['electricity_price_exp'].append(electricity_price_exp[t])
    
    #HESS
    results['P_ELY'].append(P_ELY[t])
    results['P_ElyOn'].append(P_ElyOn[t])
    results['ElyOn'].append(ElyOn[t])
    results['P_ELY_PWA'].append(P_ELY_PWA[t])
    results['mdot_H2'].append(mdot_H2[t])
    results['P_C'].append(P_C[t])
    results['E_TANK'].append(E_TANK[t])
    results['P_FC'].append(P_FC[t])
    
    # Battery
    if include_battery:
        results['E_b'].append(E_b[t])
        results['P_ch'].append(P_ch[t])
        results['P_ds'].append(P_ds[t])
    
    # Waste Heat 
    results['m_cw_ELY'].append(m_cw_ELY[t])
    results['m_cw_FC'].append(m_cw_FC[t])
    results['m_cw_HT'].append(m_cw_HT[t])
    results['m_cw_MT'].append(m_cw_MT[t])
    results['P_th_HT'].append(P_th_HT[t])
    results['P_th_MT'].append(P_th_MT[t])
    
    results['heat_zone3_35degC_demand_kWh'].append(heat_zone3_35degC_demand_kWh[t])
    results['heat_zone3_60degC_demand_kWh'].append(heat_zone3_60degC_demand_kWh[t])


# Add optimization status to results
results['status'] = m.status


print("LCOE = {:.3f} € / MWh".format(LCOE))
print("VALCOE = {:.3f} € / MWh".format(VALCOE))
print("PV Area = {:.2f} square meters".format(Area_PV))

if sum(mdot_H2) > 0:
    print("Minimal Hydrogen Tank Size = {:.2f} cubic meters".format(max(V_Tank_H2)))

#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------
# Import the plotting module

from plotting_module import (plot_power_generation, plot_component_sizes, 
                             plot_HESS_results, plot_battery_operation, 
                             plot_costs_and_prices, costs_pie_chart, plot_WHR, 
                             plot_ely_efficiency, pv_efficiency)


# Call the plotting functions as needed
# plot_power_generation(P_PV, P_imp, P_exp, df_input, nHours)
plot_component_sizes(S_PV, S_PV_max, S_ELY, S_ELY_max, S_C, S_C_max, S_FC, S_FC_max, S_TANK, S_TANK_max)
# pv_efficiency(df_pv)
plot_HESS_results(P_PV, P_ELY, P_ELY_PWA, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input)
plot_costs_and_prices(all_costs, df_input, electricity_price_imp, electricity_price_exp)
plot_WHR(results)


if include_battery:
    plot_battery_operation(P_demand, P_imp, P_ch, P_ds, E_b, bat_params, C_b_kWh, nHours)

if sum(mdot_H2) > 0:
    plot_ely_efficiency(P_ELY, P_ELY_PWA, S_ELY, nHours, x_bp_val, y_bp_val)
    
    inputPower = sorted([P_ELY[t] / S_ELY for t in range(nHours)])
    outputPower = sorted([P_ELY_PWA[t] / S_ELY for t in range(nHours)])
    eta_ELY = [(outp / inp) * 100 if inp != 0 else 0 for inp, outp in zip(inputPower, outputPower)]
    

# For plotly graphs
fig_power_generation = plot_power_generation(results, df_input, nHours)
fig_costs_pie_chart  = costs_pie_chart(all_costs)

# efficiency_fc = [P_FC_PWA[t] / u_fc[t] for t in range(nHours)] # Efficiency
# inputPower_fc = [u_fc[t] / 1000 for t in range(nHours)]
# current_fc = [i_fc[t] / 1000 for t in range(nHours)]

# # Creating the plots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# # Plotting u_fc over i_fc
# ax1.plot(inputPower_fc, current_fc, '-o', label='PWA fit')  # Convert power to kW for the plot
# ax1.set_xlabel('P [kW]')
# ax1.set_ylabel('i [A]')
# ax1.set_title('PWA approximation of the FC current')
# ax1.legend()

# # Plotting efficiency of the Fuel Cell
# ax2.plot(inputPower_fc, efficiency_fc, '-o', label='Efficiency')  # Convert power to kW for the plot
# ax2.set_xlabel('u_fc [kW]')
# ax2.set_ylabel('Efficiency')
# ax2.set_title('Efficiency of the Fuel Cell')
# ax2.legend()

# # Show the figure
# plt.tight_layout()
# plt.show()

#------------------------------------------------------------------------------
# Export results to excel 
#------------------------------------------------------------------------------

from results_export import export_optimization_results

# Define the path to the results directory
results_directory = export_path

# Export results
export_optimization_results(variable_names, results, results_directory)