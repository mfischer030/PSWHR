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
user = 'maxime'    # 'christian', 'gabriele', 'maxime', 'maxime_EMPA_WS'

from config import paths_configuration
input_path, demand_path, heat_path, function_path, export_path = paths_configuration(user)

# path to the functions directory => change in config.py if needed
sys.path.append(function_path) 

#------------------------------------------------------------------------------
# Setting up the model
#------------------------------------------------------------------------------

# Choose energy tariff
energy_tariff = "Grey" # Choose from "Green","Blue","Grey"

# Flag to include/exclude battery in the optimization model
include_battery = False

# Choose the "grid-connectivity" of the model: Grid Connectivity Factor
GCF = 10  # [%] of power which can be imported from the grid from the total energy demand

# Choose Hydrogen storage pressure in [bar] (30 bar = H2 outlet pressure from ELY => no compressor needed)
H2_storage_pressure = 30                                                     

# Capital recovery factor CRF
discountRate = 0.07                                             # [-] Annual Discount Rate; 2023_Giovanniello | 2021_Marocco

# Define number of breakpoints on PWA approx. of Electrolyser efficiency curve
N_bp = 3

# Define wether you want to consider current (1), midterm (2) or future (3) efficiencies and Unit prices for components
efficiency_level = 3
UP_level         = 3 

# Import functions and modules ------------------------------------------------

from import_data     import get_data
from efficiencies    import calculate_pv_efficiency, pwa_ELY, pwa_FC
from cost_curves     import aa_c, bb_c, xx_c, aa_e, bb_e, xx_e
from cost            import totalAnnualCost
from plotting_module import (plot_heat_demand, plot_power_generation, plot_component_sizes, 
                             plot_HESS_results, plot_battery_operation, 
                             plot_costs_and_prices, costs_pie_chart, plot_WHR, 
                             plot_efficiencies, pv_efficiency)
from results_export  import export_optimization_results

#------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------

# Import data and generate variables using the get_data function---------------
(irradiance, P_demand, T_amb, df_input, df_demand, df_heat_demand, 
 timeline_choice, season_choice) = get_data(input_path, demand_path)

# Defining heat demand variables in [Wh]
z1_35degC_kWh = df_heat_demand['Heating_Zone1_35degC_kW'].values * 1000
z1_60degC_kWh = df_heat_demand['Heating_Zone1_60degC_kW'].values * 1000
z2_35degC_kWh = df_heat_demand['Heating_Zone2_35degC_kW'].values * 1000
z2_60degC_kWh = df_heat_demand['Heating_Zone2_60degC_kW'].values * 1000
z3_35degC_kWh = df_heat_demand['Heating_Zone3_35degC_kW'].values * 1000
z3_60degC_kWh = df_heat_demand['Heating_Zone3_60degC_kW'].values * 1000

plot_heat_demand(df_heat_demand)

#------------------------------------------------------------------------------
# Input parameters
#------------------------------------------------------------------------------
# General parameters
k    = 1.4                               # Ratio cp/cv [-]
R_H2 = 4.1242 * 1000                     # Individual Gas constant H2 [J/kg*K]
M_H2 = 2.01568 / 10**-3                  # H2 molar mass in [kg/mol]
HHV  = 39.39 * 3600 * 1000               # Higher heating value of H2 in [J/kg] = 39.39 kWh/kg
W_H2_HHV = 3.545*1000                    # The vol specific chemical energy (at STP) for HHV in [Wh/m3]

# Electricity prices in [EUR/MWh]
cost_imp_el = df_input['price_Eur_MWh'].values          # hourly cost to import 
cost_exp_el = df_input['Price_DayAhed'].values          # hourly cost to export 

# Revenues from WHR in [CHF/kWh]: Information on DHN Heat prices can be found 
# here: https://www.preisueberwacher.admin.ch/pue/de/home/themen/infrastruktur/fernwaerme.html

cost_export_heatMT = 0.09     # 0.09
cost_export_heatHT = 0.1189   # 0.1189

# Time parameters
nHours = len(df_input)                        # number of hours simulated
Time   = np.arange(1, nHours + 1)             # time vector
days   = nHours / 24                          # number of days
months = df_input['MO'].unique()              # Unique months in the DataFrame
weeks  = days / 7                             # number of weeks
deltat = 3600                                 # time step (s)

kWh2J = 3600*1000                             # Convet kWh to Joule
J2kWh = 1 / (3600*1000)                       # Convert Joule to kWh 

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
epsilon         = 1e-4  

#------------------------------------------------------------------------------
# Efficiencies of components in [-]
#------------------------------------------------------------------------------
"""
PV:   from Roxanne
ELY:  2023_van_der_roest (74-79%) | IRENA: 65% | 2019_Proost Ultimate 80% | 2023_Wang: 74% | 2023_Giovanniello: 58%
C:    from Roxanne | 2020_Pan et al: 90%
TANK: 2023_Wang: 95%
FC:   2021_cigolotti: 52% 2024, 58% 2030 (approximation 2050: 70% following current trend [figure 10, c]) doi.org/10.3390/en14164963 | 2023_Wang: 50% | 2023_Giovanniello: 45% 
"""
eta_ELY_nominal = 0.65  
eta_FC_nominal  = 0.52 

# Run pwa_ELY function
x_bp_val, y_bp_val, mm_elec, qq_elec = pwa_ELY(N_bp)

eta_ELY_max = 0.7
eta_ELY_max_model = 0.555937816 # From P. Gabriellis model
scaling_eta_ELY   = eta_ELY_max / eta_ELY_max_model

y_bp_val_scaled   = y_bp_val * scaling_eta_ELY
mm_elec_scaled    = mm_elec  * scaling_eta_ELY
qq_elec_scaled    = qq_elec  * scaling_eta_ELY

if efficiency_level == 1:
    eta = {'PV': 0.21,'ELY': eta_ELY_nominal,'C': 0.7526,'TANK': 0.95,'FC': eta_FC_nominal}               # Current
elif efficiency_level == 2:
    eta = {'PV': 0.21,'ELY': 0.7,'C': 0.8,'TANK': 0.975,'FC': 0.6}                                        # Midterm
elif efficiency_level == 3: 
    eta = {'PV': 0.21,'ELY': 0.8,'C': 0.9,'TANK': 0.99,'FC': 0.7}  
    # eta = {'PV': 0.21,'ELY': 1,'C': 1,'TANK': 1,'FC': 1}                                                # Future      

#------------------------------------------------------------------------------
# Unit prices of components / capital costs in [€/W] 
#------------------------------------------------------------------------------
"""
  PV: [CHF/W]   - 2023_Tay Son Le: 0.881 | 2018_Gabrielli 300€/m2 | 2020_Pan et al: 0.9817 | 2023-30-04_Technology_characteristics_from_WP3_HSLU: 0.8
 BAT: [CHF/kWh] - 2021_Marroco
 ELY: [CHF/W]   - 2023_IEA-GlobalH2Review: PEM - 2kUSD/kW (=1841.60CHF) - reduction to 600 USD/kW (=556CHF) | 2023_ELY 1500$/kW = 1381.20 CHF/kW
   C: [CHF/W]   - 2020_Pan et al: 1228 (¥/kW) = 0.155 CHF/W | 2021_Minutillo et al 33.3 kW => 671.25 kCHF => 20.157 €/W = 19.80 CHF/W
TANK: [CHF/J]   - 1644 [EUR/kgH2] | 2018_Gabrielli: {20.7;13.6;10.9} in [€/kWh] | 2023_A Review on the Cost Analysis of H2 Gas Storage Tanks for FCV: 2020:9.18 CHF/kWh 2025: 8.25 CHF/kWh Ultimate: 7.34 CHF/kWh
  FC: [CHF/W]   - (Wang et. al 2023 - 2kUSD/kW) = 1841.43 CHF | 2018_Gabrielli {2160;1680;1320} [€/kW] | 2021_AFCTCP_Stationary_Application_Performance: Now: 2-3.5 €/W Drop: 1.2-1.75 €/W. => 1.10 CHF/W
  HP: [CHF/W]   - 2023_van der Roest: 600 €/kWth = 576 CHF/kW
 HEX: [CHF/m2]  - from Roxanne: 77.79 €/m2 (Fixed cost for heat exchanger [EUR]) | From Christian: 1782CHF/5.2m2 = 342.69 CHF/m2 => price for an actual HEX @ PSI
"""

if UP_level == 1:
    UP = {'PV': 0.8,'BAT': 0.5/3600,'ELY': 1.4,'C': 19.80, 'TANK': 9.18/(3.6*10**6),'FC': 1.84143,  'HP': 0.576, 'HEX': 342.69}        # Current
elif UP_level == 2: 
    UP = {'PV': 0.8,'BAT': 0.5/3600,'ELY': 1,'C': 14.74,'TANK': 8.25/(3.6*10**6), 'FC': 1.473145,'HP': 0.238, 'HEX': 221}              # Mid-Term
elif UP_level == 3: 
    UP = {'PV': 0.8,'BAT': 0.5/3600,'ELY': 0.556,'C': 9.82, 'TANK': 7.34/(3.6*10**6),'FC': 1.10486,'HP': 0.200, 'HEX': 100}            # Ultimate
    # UP  =  {'PV': 1,'BAT': 1,'ELY': 1,'C': 1, 'TANK': 1,'FC': 1,'HP': 1, 'HEX': 1}



# Annual maintenance cost as fraction of total cost => from Roxanne------------
maintenance = {
    'PV': 0.05,              # Annual maintenance PV - 2018_Gabrielli
    'BAT': 0.02,              # Annual cost maintenance battery - 2021_Marocco
    'ELY': 0.05,             # Annual maintenance electrolyser - 2018_Gabrielli
    'C': 0.06,               # Annual maintenance compressor - 2017_Viktorsson
    'TANK': 0.03,            # Annual maintenance storage tank - 2018_Gabrielli
    'FC': 0.08,              # Annual maintenance fuel cell - 2018_Gabrielli
    'HP': 0.015,             # Annual maintenance heat pump - 2018_Gabrielli
    'HEX': 0.02              # Annual cost maintenance HEX, frac total ann cost - 2023_vanderRoest
}

# Lifetime of the components in [years]----------------------------------------
life = {'PV': 25,    # Roxanne: 30 | 2023_Tya Son Le: 25 years
        'BAT': 10,   # From Gabriele
        'ELY': 10,   # 2023_Wang et al: 5 years | 2023_Tay Son Le: 15 years || Maxime: Before value was 20
        'C': 20,     # 2020_Pan et al: 20 years
        'TANK': 35,  # 2023_Wang et al: 20 years | 2023_Tay Son Le: 25 years | 2023-30-04_Technology_characteristics_from_WP3_HSLU: 25-50 years
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

# Maximal PV Area 
landuseCoefPV = 0.8                  # Determines how much of the disposal land can actually be used for PV power generation
Area_PV_max   = landuseCoefPV * 5000 # Maximum PV area [m2] 

# Electrolyser max and min nominal power (W)   
S_ELY_max = 1000*1000    # Maximal Electrolyzer size in [W]
S_ELY_min = 0           # Min. size ELY where problem is feasible [W] - from Rox

# Calculating the spezific work of the compresso, from Minutillo et al. 2021
# (Analyzing the levelized cost of hydrogen eq 1+2) => from Roxanne 
T_in_H2 = 65 + 273.15                                                          # H2 temperature (=T_cat=T_an) [K]  
p_out   = H2_storage_pressure                                                  # Compressor outlet pressure [bar] = H2 storage pressure 
p_in    = 30                                                                   # Compressor inlet pressure [bar]  PEM electrolyzer at Empa works at 30 bar 
L_is_C  = (k/(k-1)) * R_H2 * T_in_H2 * (((p_out/p_in)**((k-1)/k)) - 1)         # Specific work compressor [J/kg] 

# Maximal TANK energy capacity (J) => 14 days storage capacity
S_TANK_max    = E_demand_day * 14 * 3600                  # E_demand_day in [Wh]  
S_TANK_H2_max = S_TANK_max / HHV                          # equivalent in kg_H2

# Maximal FC size as the maximal power demand divided by eff in [W]
#S_FC_max = P_peak_max / eta["FC"]
S_FC_max = 1000*1000  
"""
Range PEMFC = [10W;1MW] from 2021_cigolotti Comprehensive Review on Fuel Cell 
Technology for Stationary Applications
"""
# Waste Heat Recovery ---------------------------------------------------------

T_in  = 57                     # Inlet temperature cooling water to HEX in [°C]
T_out = 62                     # Temperature cooling water to applications [°C]
T_HEX = 64                     # Outlet temperature cooling water from HEX [°C]
c_p   = 4186                   # Specific heat capacity cooling water [J/kgK]

P_th_ELY_max = (1 - eta['ELY']) * S_ELY_max # max heat recovered from ELY in [W]
P_th_FC_max  = (1 - eta['FC'])  * S_FC_max  # max heat recovered from FC in [W]
P_th_max     = P_th_ELY_max + P_th_FC_max

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
COP_carnot = T_HT_out / (T_HT_out - T_out)             # Maximum COP HP
COP        = 0.5 * COP_carnot                          # Real COP, as in Tiktak

#------------------------------------------------------------------------------

# Prompt the user to choose "Grid" or "Off-Grid" scenario
scenario_choice = input("Enter 'grid' for grid-connected scenario or 'off-grid' for an off-grid scenario: ").lower()
# Check the user's choice and set relevant parameters accordingly
if scenario_choice == 'grid':
    P_imp_ub = 10000000                 # Upper bound for P_imp in grid-connected scenario, Trafo limit is 1MW
    P_exp_ub = 1000000                  # Upper bound for P_exp in grid-connected scenario, Trafo limit is 1MW
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
Area_PV  = m.addVar(lb=0,         ub=Area_PV_max, name='Area_PV')              # Maximum PV Area in [m2]
S_ELY    = m.addVar(lb=S_ELY_min, ub=S_ELY_max,   name='S_ELY')                # Size of the Electrolyzer in [W]
S_C      = m.addVar(lb=0,         ub=100000,      name='S_C')                  # Size of the Compressor in [W]
S_TANK   = m.addVar(lb=0,         ub=S_TANK_max,  name='S_TANK')               # Size of the Hydrogen Storage Tank in [J]
S_FC     = m.addVar(lb=0,         ub=S_FC_max,    name='S_FC')                 # Size of the Fuel CEll in [W]
S_HP     = m.addVar(lb=0, ub=P_th_max, name="S_HP")                            # Size of the Heat Pump in [W]
S_HEX    = m.addVar(lb=0, ub=S_HEX_max, name="S_HEX")                          # Size of the Heat Exchanger in [m2]

# Operation Variables 
P_ELY     = m.addVars(nHours, lb=0, ub=S_ELY_max, name='P_ELY')                # Electrolyser input power in [W]
# P_ELY_out = m.addVars(nHours, lb=0, ub=S_ELY_max, name='P_ELY')                # Output from the Electrolyser => P_ELY[t] * eta[ELY][t]
# P_ElyOn   = m.addVars(nHours, lb=0, ub=S_ELY_max, name="P_ElyOn")
# ElyOn     = m.addVars(nHours, lb=0, ub=1,vtype=GRB.INTEGER, name="ElyOn")      # Binary Variable for ON-OFF condition of the Electrolyser

E_TANK    = m.addVars(nHours, lb=0, ub=S_TANK_max, name='E_TANK')              # Hydrogen Energy stored in Tank in [J]

P_FC       = m.addVars(nHours, lb=0, ub=S_FC_max, name='P_FC')                 # Fuel Cell Output Power in [W]
# i_FC       = m.addVars(nHours, name="i_FC",       vtype=GRB.CONTINUOUS)        # Fuel Cell current in [A]
# Vdot_FC_H2 = m.addVars(nHours, name="Vdot_FC_H2", vtype=GRB.CONTINUOUS)        # FC inlet volume flow 
# P_FC_in    = m.addVars(nHours, name="P_FC_in",    vtype=GRB.CONTINUOUS)        # Fuel Cell Input Power in [W]
    
P_imp     = m.addVars(nHours, lb=0, ub=P_imp_ub,   name="P_imp")               # Imported electricity from the Grid in [W]
P_max_imp = m.addVars(months, lb=0, ub=P_imp_ub,   name="P_max_imp")           # Variable introduce to define the maximum grid import during a month (needed to calculate the Benutzungsdauer which determines the grid use tariff)
P_exp     = m.addVars(nHours, lb=0, ub=P_exp_ub,   name="P_exp")               # Exported electricity to the Grid in [W] # Check upper-bound; 

# Waste Heat recovery varibles
m_cw_ELY = m.addVars(nHours, lb=0, ub=m_cw_ELY_max, name="m_cw")               # Cooling water ELY mass flow in [?]
m_cw_FC  = m.addVars(nHours, lb=0, ub=m_cw_FC_max, name="m_cw")                # Cooling water ELY mass flow in [?]

m_cw_HT  = m.addVars(nHours, lb=0, ub=m_cw_max, name="m_cw_HT")                # Hight Temp water mass flow ELY + FC
m_cw_MT  = m.addVars(nHours, lb=0, ub=m_cw_max, name="m_cw_LT")                # Medium Temp. water mass flow ELY + FC
P_th_HT  = m.addVars(nHours, lb=0, ub=P_th_max, name="P_th_HT")  
P_th_MT  = m.addVars(nHours, lb=0, ub=P_th_max, name="P_th_LT")


if include_battery:
    P_ch = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="P_ch")             # Power charged to the battery
    P_ds = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="P_ds")             # Power discharged from the battery
    E_b  = m.addVars(nHours,lb=0, ub=bat_params['C_b_max'], name="E_b")              # Energy in the battery at each time step
    C_b  = m.addVar(lb=bat_params['C_b_min'], ub=bat_params['C_b_max'], name="C_b")  # Battery Capacity
    

grid_usage = m.addVars(nHours, vtype=GRB.BINARY, name="grid_usage")            # Binary variable for grid usage indicator, 1 if used, i.e. P_imp>0, 0 otherwise
h_usage    = m.addVar(lb=0, ub=nHours, name="h_usage")                         # Total hours of grid usage
high_usage = m.addVar(vtype=GRB.BINARY, name="high_usage")                     # Binary variable indicating whether grid usage is above the threshold, 1 if above, 0 if below

# Variables
BD = m.addVar(vtype=GRB.CONTINUOUS, name="BD")


#------------------------------------------------------------------------------
# Additional calculations
#------------------------------------------------------------------------------

# PV calculation---------------------------------------------------------------

# Calling the function which calculates the relative efficiency
eta_cell = calculate_pv_efficiency(irradiance, T_amb)

# Generate dataframe for plot
df_pv = pd.DataFrame({'irradiance': irradiance,'T_amb': T_amb,'eta_cell': eta_cell})

P_PV = [irradiance[t] * eta_cell[t] * Area_PV  for t in range(nHours)]          # Hourly PV power generation [W]
    
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
mdot_H2 = [P_ELY[t] * eta['ELY']  / HHV for t in range(nHours)] 
#mdot_H2   = [P_ELY_PWA[t] / HHV for t in range(nHours)] 

# Size (S_C) and operating power (P_C) of the compressor in [W]
m.addConstr(S_C == mdot_H2_nom * L_is_C / eta["C"], name='HESS Balance')

# S_C     = mdot_H2_nom * L_is_C / eta["C"]  
S_C_max = mdot_H2_max * L_is_C / eta["C"] 
P_C     = [mdot_H2[t] * L_is_C / eta['C'] for t in range(nHours)]

#------------------------------------------------------------------------------
# Constraints
#------------------------------------------------------------------------------

for t in range(1, nHours):
    # Energy Balance for HESS in [J] considering power a timestep t in [W]
    m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - (P_FC[t] * deltat) / eta['FC'], name='HESS Balance') # Before 
    # m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - (P_FC_in[t] * deltat), name='HESS Balance') 
    
    # Constraint for the ELY output power not to exceed the storage capacity
    # m.addConstr(P_ELY[t]  <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY') 
    m.addConstr(P_ELY[t] * eta["ELY"] <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY')       
    
    # Constraint for the FC to only use hydrogen available in the storage
    # m.addConstr(P_FC[t] <= E_TANK[t-1] / deltat, name='FC')
    m.addConstr(P_FC[t] / eta['FC'] <= E_TANK[t-1] / deltat, name='FC')  
                  
for t in range(nHours):
    m.addConstr((P_ELY[t]  <= S_ELY),  name="upper_Size_Constraint_ELY")
    m.addConstr((E_TANK[t] <= S_TANK), name= "upper_Size_Constraint_TANK")
    m.addConstr((P_FC[t] / eta['FC'] <= S_FC ), name= "upper_Size_Constraint_FC")
    

# Initializing FC power
m.addConstr(P_FC[0] <= E_TANK[0] / deltat, name= "InitialFC")                                    

# constraint for H2 storage equal at final and last time step (periodicity)
m.addConstr(E_TANK[0] == E_TANK[nHours-1], name='Periodicity_HESS') # 08.02: added name to cosntraint

# Overall energy balance: left => consumers | right => generators
if include_battery:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t] + P_ch[t] + P_exp[t] 
                  + P_demand[t] <= P_PV[t] + P_ds[t] + P_imp[t] + P_FC[t]  
                  for t in range(nHours)), name='EnergyBalance') 
else:
    m.addConstrs((P_ELY[t] + (P_th_HT[t])/COP + P_C[t] + P_exp[t] + 
                  P_demand[t] <= P_PV[t] + P_imp[t] + P_FC[t] 
                  for t in range(nHours)), name='EnergyBalance')

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
        m.addConstr(E_b[t] == E_b[t-1] * bat_params['eff_sd'] + P_ch[t] * 
                    bat_params['eff_ch'] * deltat 
                    - ((P_ds[t] * deltat) / bat_params['eff_disch']), 
                    name='Battery energy balance')
    
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
    m.addConstr(m_cw_ELY[t] ==  ( (1 - eta['ELY']) * P_ELY[t])    / (c_p * (T_HEX - T_in)), "PEM_outlet")
    m.addConstr(m_cw_FC[t]  ==  ( ((1/eta['FC']) - 1) *  P_FC[t])  / (c_p * (T_HEX - T_in)), "FC_outlet")

    
    # Cooling flow requirements
    m.addConstr(m_cw_ELY[t] + m_cw_FC[t] >= m_cw_HT[t] + m_cw_MT[t], "massflowBalanceCoolingWater")
    
    # Low- and High-temperature heat output
    m.addConstr(P_th_MT[t] == m_cw_MT[t] * c_p * (T_out - T_in), "LTheat")
    m.addConstr(P_th_HT[t] == m_cw_HT[t] * c_p * (T_out - T_in), "HTheat")
    
    # Heat exchanger size constraint for medium-temperature side
    m.addConstr(P_th_MT[t] / (U_HEX * T_log) <= S_HEX, "HEXsize")
    
    # Heat pump size constraint for high-temperature side
    m.addConstr(P_th_HT[t] <= S_HP, "HPsize")
    
    # Heat demand constraint
    m.addConstr(P_th_MT[t] <= z1_35degC_kWh[t] + z2_35degC_kWh[t] + z3_35degC_kWh[t], "Heat demand satisfaction") # Multiply by 1000 to conver kW to W
    m.addConstr(P_th_HT[t] <= z1_60degC_kWh[t] + z2_60degC_kWh[t] + z3_60degC_kWh[t], "Heat demand satisfaction")

# BIG-M constraint for Benutzungsdauer = grid usage time-----------------------
# # we use the big-M approach to define the binary variable grid_usage
# for i in range(nHours):
#     m.addConstr(P_imp[i] <= M_hours * grid_usage[i])                           # here we constraint grid_usage to take value 1 if P_imp is larger than zero. If P_imp is positive, this constraint is satisfied only if grid_usage is 1
#     m.addConstr(P_imp[i] >= epsilon - M_hours * (1 - grid_usage[i]))           # here we constraint grid_usage to 0 if P_imp is zero. 
# # Here I choose M_hours in the order of the peak demand

# # here we add the constraints needed for the grid usage tarif selection
# m.addConstr(h_usage == gp.quicksum(grid_usage[i] for i in range(nHours)))      # not sure this has to be a constraint, could be a normal calculation maybe 

# # defining the value of the binary variable high-usage based on the hours of operation, h_usage
# m.addConstr(h_usage - threshold_hours[timeline_choice] <= M_threshold * high_usage)# if h_usage is larger than threshold, the inequality is respected only if high_usage takes value 1. 
# # here I choose the M_threshold in the order of the hours in the year

# imported power and average of imported peak
total_imported_power = gp.quicksum(P_imp[i] for i in range(nHours))
avg_monthly_max      = gp.quicksum(P_max_imp) / len(months)


# # BD constraints using binary variable and bigM
# m.addConstr(total_imported_power - threshold_hours[timeline_choice] * avg_monthly_max >= -M_hours * (1 - high_usage), "activate_high_usage")
# m.addConstr(total_imported_power - threshold_hours[timeline_choice] * avg_monthly_max <= M_hours  * high_usage, "deactivate_high_usage")

# Constraint: Calculate BD
m.addConstr(BD *  avg_monthly_max == total_imported_power, "Calculate_BD")

# # Big M method constraints
M = 10000
epsilon=0.001

# m.addConstr(BD - 3500 >= -M * (1 - high_usage), "BD_Less_Than_3500")
# m.addConstr(BD - 3500 <= M * high_usage, "BD_More_Than_3500")

# PWA approximation of the Electrolyser----------------------------------------

# for t in range(nHours):
    
#     m.addGenConstrPWL(P_ELY[t], P_ELY_PWA[t], x_bp_val*100*1000, y_bp_val_scaled*100*1000)
    
    # # Constraint: P_e should be less than or equal to P_ElyOn
    # m.addConstr(P_ELY[t] <= P_ElyOn[t], "MaxPowerEly")

    # # Constraint: P_e should be at least 20% of P_ElyOn
    # m.addConstr(P_ELY[t] >= 0.2 * P_ElyOn[t], "MinPowerEly")

    # # Constraint: P_ElyOn not to exceed the maximum power when the electrolyser is on (S_e_max) 
    # m.addConstr(P_ElyOn[t] <= S_ELY_max * ElyOn[t], "PowerEly1")

    # # Constraint: P_ElyOn should not exceed S_e
    # m.addConstr(P_ElyOn[t] <= S_ELY, "PowerEly3")

    # # Constraint: When the electrolyser is off (ElyOn = 0), P_ElyOn should be at least S_e - S_e_max. When ElyOn = 1, this constraint has no effect because S_e - S_e_max * (1 - ElyOn) equals S_e - S_e_max * 0 = S_e
    # m.addConstr(P_ElyOn[t] >= S_ELY - S_ELY_max * (1 - ElyOn[t]), f"PowerEly4_{i}")
    
    # # Piecewise linear approximation constraints
    # if N_bp == 1:
    #     m.addConstr(P_ELY_PWA[t] == eta_ELY_nominal * P_ELY[t], "Efficiency Electrolyser N_bp=1")
    # elif N_bp == 2:
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled * P_ELY[t] + qq_elec_scaled * P_ElyOn[t],  "PWA_Q1")
    # elif N_bp==3:
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled[0] * P_ELY[t] + qq_elec_scaled[0] * P_ElyOn[t], "PWA_Q1")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled[1] * P_ELY[t] + qq_elec_scaled[1] * P_ElyOn[t], "PWA_Q2")
    # elif N_bp==4:
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled[0] * P_ELY[t] + qq_elec_scaled[0] * P_ElyOn[t], "PWA_Q1")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled[1] * P_ELY[t] + qq_elec_scaled[1] * P_ElyOn[t], "PWA_Q2")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec_scaled[2] * P_ELY[t] + qq_elec_scaled[2] * P_ElyOn[t], "PWA_Q3")
    # elif N_bp==8:
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[0]  * P_ELY[t] + qq_elec[0] * P_ElyOn[t],  "PWA_Q1")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[1]  * P_ELY[t] + qq_elec[1] * P_ElyOn[t], "PWA_Q2")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[2]  * P_ELY[t] + qq_elec[2] * P_ElyOn[t], "PWA_Q3")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[3]  * P_ELY[t] + qq_elec[3] * P_ElyOn[t],  "PWA_Q4")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[4]  * P_ELY[t] + qq_elec[4] * P_ElyOn[t],  "PWA_Q5")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[5]  * P_ELY[t] + qq_elec[5] * P_ElyOn[t],  "PWA_Q6")
    #     m.addConstr(P_ELY_PWA[t] <= mm_elec[6]  * P_ELY[t] + qq_elec[6] * P_ElyOn[t],  "PWA_Q7")
        
# PWA approximation of the Fuel Cell-------------------------------------------

# # Run the pwa_FC function to get FC parameters
# FC_param_struct = pwa_FC()

# # Adding PWA constraints for each hour
# breakpoint_x = FC_param_struct['FC_p_breakpoint']
# breakpoint_y = FC_param_struct['FC_i_breakpoint']
# slope_low    = FC_param_struct['FC_slope_low_current']
# slope_high   = FC_param_struct['FC_slope_high_current']

# for t in range(nHours):
#     # Piecewise constraint for i_FC[t]
#     m.addGenConstrPWL(P_FC[t], i_FC[t], [0, breakpoint_x, FC_param_struct['P_FC_max']],
#                       [0, breakpoint_y, slope_high * (FC_param_struct['P_FC_max'] - breakpoint_x) + breakpoint_y])

#     # Constraint to convert i_FC[t] to Vdot_FC_H2[t]
#     m.addConstr(Vdot_FC_H2[t] == FC_param_struct['FC_i_to_VH2'] * i_FC[t], name=f"VH2_conversion_{t}")

#     # Constraint for P_FC_in[t] calculation
#     m.addConstr(P_FC_in[t] == W_H2_HHV * Vdot_FC_H2[t], name=f"power_input_conversion_{t}")

# Implement cost curves: compressor and electrolyser --------------------------

# Define the breakpoints for the piecewise approximation of the cost curves
x_bp_ELY = [0, 250000, 400000, 550000]
y_bp_ELY = [0, 347388, 495479, 632884]

x_bp_C = [0, 5000, 12500, 20000]
y_bp_C = [0, 112682, 192791, 253934]

# Define binary variables for selecting the piecewise segments for ELY and C
y_ELY = m.addVars(len(x_bp_ELY) + 1, vtype=GRB.BINARY, name="y_ELY")
y_C   = m.addVars(len(x_bp_C) + 1, vtype=GRB.BINARY, name="y_C")

# Ensure only one y variable for ELY and C can be 1 at the same time
m.addConstr(sum(y_ELY[i] for i in range(len(y_ELY))) == 1, "UniqueSegment_ELY")
m.addConstr(sum(y_C[i]   for i in range(len(y_C)))   == 1, "UniqueSegment_C")

# Add the piecewise linear cost functions for ELY using the breakpoints
cost_ELY = m.addVar(vtype=GRB.CONTINUOUS, name="cost_ELY")

m.addGenConstrPWL(S_ELY, cost_ELY, x_bp_ELY, y_bp_ELY, "PWL_ELY")

# # Add the piecewise linear cost functions using the breakpoints for C
cost_C = m.addVar(vtype=GRB.CONTINUOUS, name="cost_C")
m.addGenConstrPWL(S_C, cost_C, x_bp_C, y_bp_C, "PWL_C")

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
     
[cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, 
  cost_maint, cost_WHR, all_costs] = totalAnnualCost(
                        system_sizes, energy_tariff, discountRate,
                        UP, maintenance, life, 
                        P_imp, P_max_imp, P_exp, P_th_MT, P_th_HT,
                        cost_imp_el, cost_exp_el, cost_export_heatMT, cost_export_heatHT,
                        m, high_usage,
                        df_input, nHours, timeline_choice,
                        cost_ELY, cost_C,
                        UP_level
                        )

# Total annual costs in [€/y]--------------------------------------------------

# Defining a factor to adjust overall cost to choosen timeline
if timeline_choice == 'week':
    multiplier = 52                               # Assuming 52 weeks in a year
elif timeline_choice == 'month':
    multiplier = 12                               # 12 months in a year
else:
    multiplier = 1                                # Default to yearly calculation with no multiplication needed

cost = (cost_inst + cost_maint)/multiplier + cost_op

all_costs['TAC'] = cost

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

# For Gurobi tupledict object--------------------------------------------------
P_ELY      = [P_ELY[t].X   for t in range(nHours)]
# P_ElyOn    = [P_ElyOn[t].X for t in range(nHours)]
# P_ELY_PWA  = [P_ELY_PWA[t].X for t in range(nHours)]
# ElyOn      = [ElyOn[t].X   for t in range(nHours)]

P_FC       = [P_FC[t].X for t in range(nHours)]
# i_FC       = [i_FC[t].X for t in range(nHours)]
# P_FC_in    = [P_FC_in[t].X for t in range(nHours)]
# Vdot_FC_H2 = [Vdot_FC_H2[t].X for t in range(nHours)]

E_TANK  = [E_TANK[t].X  for t in range(nHours)]
P_imp   = [P_imp[t].X   for t in range(nHours)]
P_exp   = [P_exp[t].X   for t in range(nHours)]
P_max_imp = [P_max_imp[month].X for month in months]

if include_battery:
    E_b     = [E_b[t].X  for t in range(nHours)]
    P_ch    = [P_ch[t].X for t in range(nHours)]
    P_ds    = [P_ds[t].X for t in range(nHours)]
    C_b     = C_b.X
    S_BAT   = S_BAT.X

m_cw_ELY = [m_cw_ELY[t].X  for t in range(nHours)]
m_cw_FC  = [m_cw_FC[t].X   for t in range(nHours)]
m_cw_HT  = [m_cw_HT[t].X   for t in range(nHours)]
m_cw_MT  = [m_cw_MT[t].X   for t in range(nHours)]
P_th_HT  = [P_th_HT[t].X   for t in range(nHours)]
P_th_MT  = [P_th_MT[t].X   for t in range(nHours)]

# For lists with Gurobi LinExpr------------------------------------------------
P_PV    = [P_PV[t].getValue()    for t in range(nHours)]
P_C     = [P_C[t].getValue()     for t in range(nHours)]
mdot_H2 = [mdot_H2[t].getValue() for t in range(nHours)]

# For Gurobi var object--------------------------------------------------------
Area_PV    = Area_PV.X
S_ELY      = S_ELY.X
S_C        = S_C.X
S_TANK     = S_TANK.X
S_FC       = S_FC.X
S_HEX      = S_HEX.X
S_HP       = S_HP.X
high_usage = high_usage.X
h_usage    = h_usage.X

cost_C   = cost_C.X
cost_ELY = cost_ELY.X
BD      = BD.X

# For Gurobi LinExpr-----------------------------------------------------------
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
    """
    Recursively retrieve values from a nested dictionary of Gurobi Var, LinExpr objects, and floats.

    Parameters:
    - costs: Dictionary containing costs as values, which may be floats, LinExpr objects, or Var objects.

    Returns:
    - Dictionary with the same structure where Gurobi objects are replaced by their numerical values.
    """
    retrieved_costs = {}
    for key, value in costs.items():
        if isinstance(value, dict):
            # Recurse into the sub-dictionary
            retrieved_costs[key] = retrieve_values(value)
        elif hasattr(value, 'X'):
            # If the object is a Gurobi Var, use the .X attribute to get the value
            retrieved_costs[key] = value.X
        elif hasattr(value, 'getValue'):
            # Check if the object has the 'getValue' method typical of Gurobi LinExpr objects
            retrieved_costs[key] = value.getValue()
        else:
            # Handle floats and other types that do not need special methods to retrieve the value
            retrieved_costs[key] = value
    return retrieved_costs

# Assuming costs_with_annuity is defined and some entries might be Gurobi LinExpr objects or floats
all_costs    = retrieve_values(all_costs)
system_sizes = retrieve_values(system_sizes)

#------------------------------------------------------------------------------
# Post processing
#------------------------------------------------------------------------------
#Calculate the Benutzungsdauer
Benutzungsdauer = sum(P_imp)/np.mean(P_max_imp)


S_PV_max = 1000*eta['PV']*Area_PV_max
S_C_max = max(mdot_H2) * L_is_C / (eta["C"] * deltat)

# Convert the battery capacity from Joules to Kilowatthours
if include_battery:
    C_b_kWh = C_b * J2kWh 

# Calculate the electricity price from BKW ------------------------------------
energy_price = all_costs["operational cost"]["electricity prices"]['Electricity prices [Rp./kWh]']
grid_fee     = all_costs["operational cost"]["electricity prices"]['Grid use cost [Rp./kWh]']
taxes        = all_costs["operational cost"]["electricity prices"]['Taxes & Levies [Rp/kWh]']
export_price = all_costs["operational cost"]["electricity prices"]['Elecricity export price [Rp.kWh]']

# Import and Export prices in CHF/MWh
electricity_price_imp = [(ip + gf + tx) * 10 for ip, gf, tx in zip(energy_price, grid_fee, taxes)]
electricity_price_exp = [ep             * 10 for ep in export_price]
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

LCOE   = cost / (sum(P_imp + P_PV)/10**6)  # Calculate LCOE # Look into the literature how
VALCOE = cost / (sum(P_demand)/10**6)

#------------------------------------------------------------------------------
# Display results
#------------------------------------------------------------------------------
simulation_report = {
    "Choosen Timeline": timeline_choice,
    "Chosen energy tariff": energy_tariff,
    "Include battery in the simulation": include_battery,
    "Grid Connectivity Factor (GCF) (%)": GCF,
    "Components efficiencies": efficiency_level,
    "Components unit prices": UP_level,
    "Maximal PV Area (m2)": Area_PV_max,
    "Maximal ELY Capacity (kW)": S_ELY_max/1000,
    "Maximal FC Capacity (kW)": S_FC_max/1000,
    "H2 Storage Pressure (bar)": p_out,
    "LCOE in CHF/MWh": f"{LCOE:.2f}",
    "Benutzungsdauer": f"{Benutzungsdauer:.2f}"
}

# Convert dictionary to DataFrame for better visualization
df_simulation_report = pd.DataFrame(list(simulation_report.items()), columns=['Parameter', 'Value'])
# print(df_simulation_report)

# Create a figure and an axes to display the table with adjusted settings
fig, ax = plt.subplots(figsize=(12, 4))  # Increased figure size for better content fit
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=df_simulation_report.values, colLabels=df_simulation_report.columns, loc='center', cellLoc='right')
table.auto_set_font_size(False)
table.set_fontsize(16)  # Slightly reduced font size
table.scale(1, 2)  # Adjusted cell scaling for better visibility and to prevent overlap
plt.title('Simulation Report')
plt.show()



# Generate a dictionnary with the results
if include_battery:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
        'P_exp', 'electricity_price_exp', 
        'P_ELY', 'mdot_H2', 'P_C', 'E_TANK','P_FC', 
        'E_b', 'P_ch', 'P_ds', 
        'm_cw_ELY', 'm_cw_FC', 'm_cw_HT','m_cw_MT', 'P_th_HT', 'P_th_MT', 
        'z1_35degC_kWh', 'z1_60degC_kWh', 'z2_35degC_kWh', 'z2_60degC_kWh',
        'z3_35degC_kWh', 'z3_60degC_kWh'
        ]

else:
    variable_names = [
        'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
        'P_exp', 'electricity_price_exp', 
        'P_ELY', 'mdot_H2', 'P_C', 'E_TANK', 'P_FC', 
        'm_cw_ELY', 'm_cw_FC', 'm_cw_HT', 'm_cw_MT', 'P_th_HT', 'P_th_MT', 
        'z1_35degC_kWh', 'z1_60degC_kWh', 'z2_35degC_kWh', 'z2_60degC_kWh',
        'z3_35degC_kWh', 'z3_60degC_kWh'
        ]

# else:
#     variable_names = [
#         'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
#         'P_exp', 'electricity_price_exp', 
#         'P_ELY', 'mdot_H2', 'P_C','E_TANK', 'P_FC_in', 'P_FC', 'i_FC', 'Vdot_FC_H2',
#         'm_cw_ELY', 'm_cw_FC', 'm_cw_HT', 'm_cw_MT', 'P_th_HT', 'P_th_MT',
#         'heat_zone3_35degC_demand_kWh', 'heat_zone3_60degC_demand_kWh'
#         ]
    
# else:
#     variable_names = [
#         'irradiance', 'P_demand', 'P_PV', 'P_imp', 'electricity_price_imp',
#         'P_exp', 'electricity_price_exp', 
#         'P_ELY', 'P_ElyOn', 'ElyOn', 'P_ELY_PWA', 'mdot_H2', 'P_C', 
#         'E_TANK', 'P_FC_in', 'P_FC', 'i_FC', 'Vdot_FC_H2',
#         'm_cw_ELY', 'm_cw_FC', 'm_cw_HT', 'm_cw_MT', 'P_th_HT', 'P_th_MT',
#         'heat_zone3_35degC_demand_kWh', 'heat_zone3_60degC_demand_kWh'
#         ]

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
    # results['P_ElyOn'].append(P_ElyOn[t])
    # results['ElyOn'].append(ElyOn[t])
    # results['P_ELY_PWA'].append(P_ELY_PWA[t])
    results['mdot_H2'].append(mdot_H2[t])
    results['P_C'].append(P_C[t])
    results['E_TANK'].append(E_TANK[t])
   
    # results['P_FC_in'].append(P_FC_in[t])
    results['P_FC'].append(P_FC[t])
    # results['i_FC'].append(i_FC[t])
    # results['Vdot_FC_H2'].append(Vdot_FC_H2[t])
    
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
    
    results['z1_35degC_kWh'].append(z1_35degC_kWh[t])
    results['z1_60degC_kWh'].append(z1_60degC_kWh[t])
    results['z2_35degC_kWh'].append(z1_35degC_kWh[t])
    results['z2_60degC_kWh'].append(z1_60degC_kWh[t])
    results['z3_35degC_kWh'].append(z1_35degC_kWh[t])
    results['z3_60degC_kWh'].append(z1_60degC_kWh[t])


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
# Call the plotting functions as needed

plot_power_generation(results, df_input, nHours)
plot_component_sizes(S_PV, S_PV_max, S_ELY, S_ELY_max, S_C, S_C_max, S_FC, S_FC_max, S_TANK, S_TANK_max)
# pv_efficiency(df_pv)

plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input) 
# plot_HESS_results(P_PV, P_ELY, P_ELY_PWA, S_ELY, S_ELY_max, P_FC, P_FC_in, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input) # if PWA

plot_costs_and_prices(all_costs, df_input, electricity_price_imp, electricity_price_exp)
plot_WHR(results)

if include_battery:
    plot_battery_operation(P_demand, P_imp, P_ch, P_ds, E_b, bat_params, C_b_kWh, nHours)

# if sum(mdot_H2) > 0:                                                         # if PWA
#     plot_efficiencies(P_ELY, P_ELY_PWA, S_ELY, nHours, x_bp_val, y_bp_val_scaled, 
#                       Vdot_FC_H2, i_FC, P_FC_in, P_FC, S_FC)
    

# For plotly graphs
# fig_power_generation = plot_power_generation(results, df_input, nHours)
# fig_costs_pie_chart  = costs_pie_chart(all_costs)

#------------------------------------------------------------------------------
# Export results to excel 
#------------------------------------------------------------------------------

# Define the path to the results directory
results_directory = export_path

# Export results
export_optimization_results(variable_names, results, results_directory, all_costs, system_sizes)