# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 18:26:35 2023

@author: fism
"""
import pandas as pd
import hvplot.pandas
import numpy as np
import matplotlib.pyplot as plt
from gurobipy import Model, GRB
import sys # 
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
# Import data
#------------------------------------------------------------------------------
# path to the module 'import_data' => change path if needed
sys.path.append(r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\functions') 

# Import data and generate variables using the get_data function
from import_data import get_data
df_input, df_demand, irradiance, P_demand = get_data() 

#------------------------------------------------------------------------------
# Input parameters
#------------------------------------------------------------------------------
# General parameters
k    = 1.4                               # Ratio cp/cv [-]
R_H2 = 4.1242                            # Individual Gas constant H2 [kJ/kg*K]
HHV  = 39.39*3600*10**3                  # Higher heating value of H2 in [J/kg]

project_lifetime     = 25        # in [y] not used for the moment...
"""Gabriele: do you have a reference for this value by chance?
   Maxime: 2023_Wang et al - 20 years
   Define the project life dependent on the smallest component lifetime?
"""
annual_interest_rate = 0.04      # Discount rate, as encouraged by EU, from Rox


# Electricity prices in [EUR/MWh]
cost_imp_el = df_input['price_Eur_MWh'].values   # hourly cost to import Maxime: values changed
cost_exp_el = df_input['Price_DayAhed'].values       # hourly cost to export 

# Time parameters
nHours = len(df_input)                             # number of hours simulated
Time   = np.arange(1, nHours + 1)                  # time vector
days   = nHours / 24                               # number of days
weeks  = days / 7                                  # number of weeks
deltat = 3600                                      # time step (s)

# Efficiencies of components in [-] 
# eta = {'PV': 0.21,     # from Roxanne
#        'ELY': 0.6,     # 
#        'C': 0.7526,    # from Roxanne
#        'TANK': 0.95,   # 2023_Wang et al
#        'FC': 0.5       #
# }     

eta = {'PV': 0.21,     # from Roxanne
       'ELY': 1,     # 
       'C': 1,         # from Roxanne
       'TANK': 1,      # 2023_Wang et al
       'FC': 1         #
}     
"""
2021_cigolotti Comprehensive Review on Fuel Cell Technology for Stationary 
Applications:
    Electric Efficiency PEMFC: [38,38,37,40] in [%]
"""
# Lifetime of the components in [years]
life = {'PV': 25,    # Roxanne: 30 | 2023_Tya Son Le: 25 years
        'ELY': 10,   # 2023_Wan et al: 5 years | 2023_Tya Son Le: 15 years || Maxime: Before value was 20
        'C': 20,     # 
        'TANK': 35,  # 2023_Wan et al: 20 years | 2023_Tya Son Le: 25 years
        'FC': 5      # 2023_Wan et al: 5 years | 2023_Tya Son Le: 5 years || Maxime: Before value was 20
}
"""
#Gabriele: do you have references for these values? Can you please check 
the lifetime values for ely and FC
"""
# Unit prices of components / capital costs in [€/W] => from Roxanne
# UP = {
#     'PV': 0.8,        # Unit price PV in [€/W] | 2023_Tya Son Le: 881 USD/year | 2018_Gabrielli 300€/m2 
#     'ELY': 1.7,       # [€/W] (2023_IEA-GlobalH2Review: PEM - 2kUSD/kW)
#     'C': 0.5,         # CHECK FOR DATA
#     'TANK': 1644/HHV, # Unit price hydrogen storage [EUR/J], 1644 [EUR/kgH2] | 2018_Gabrielli: {20.7;13.6;10.9} in [€/kWh]
#     'FC': 1.680,      # [€/W] (Wang et. al 2023 - 2kUSD/kW) | 2018_Gabrielli {2160;1680;1320} [€/kW]
#     'HP': 0.576,      # Unit price heat pump
# }

UP = {
    'PV': 0.8,           # Unit price PV in [€/W] | 2023_Tya Son Le: 881 USD/year | 2018_Gabrielli 300€/m2 
    'ELY': 1.7/10,       # [€/W] (2023_IEA-GlobalH2Review: PEM - 2kUSD/kW)
    'C': 0.5/10,         # CHECK FOR DATA
    'TANK': 1644/(10*HHV), # Unit price hydrogen storage [EUR/J], 1644 [EUR/kgH2] | 2018_Gabrielli: {20.7;13.6;10.9} in [€/kWh]
    'FC': 1.680/10,      # [€/W] (Wang et. al 2023 - 2kUSD/kW) | 2018_Gabrielli {2160;1680;1320} [€/kW]
    'HP': 0.576,      # Unit price heat pump
}
# Annual maintenance cost as fraction of total cost => from Roxanne
maintenance = {
    'PV': 0.0158,       # Annual maintenance PV
    'ELY': 0.02,        # Annual maintenance electrolyser
    'C': 0.08,          # Annual maintenance compressor
    'TANK': 0.03,       # Annual maintenance storage tank
    'FC': 0.02,         # Annual maintenance fuel cell
    'HP': 0.015,        # Annual maintenance heat pump
}
# Daily energy demand in [Wh/day] - P_demand is in [W]
E_demand_day = sum(P_demand)/(4*days) # Gabriele: this is [Wh], correct | Maxime: sum(P_demand) is Wh divided by number of days => Wh/day

# Maximal daily power demand in [W]
P_peak_max = np.max(P_demand)

# Defining maximal SIZES of the components-------------------------------------
# maximum PV size needs to be defined to define the upper bound of the design variable.
# Area_PV_max  = P_peak_max / (eta['PV'] * np.mean(irradiance)) # Maximum PV area [m2]
Area_PV_max  = 5000 # Maximum PV area [m2]

# Electrolyser maximal nominal power (W) very large (not realistic)  
S_ELY_max = 500*1000                         # Maximal Electrolyzer size in [W]

# Calculating the spezific work of the compresso, from Minutillo et al. 2021
# (Analyzing the levelized cost of hydrogen eq 1+2) => from Roxanne 
T_in_H2 = 70 + 273.15                       # H2 temperature (=T_cat=T_an) [K]  # Gabriele: this could probably be 65 degC, see for example: https://www.sciencedirect.com/science/article/pii/S0360319923015410
p_out   = 350                               # Compressor outlet pressure [bar]  # Gabriele: is your p_out at 820 bar? Please see in the litearture what value is set here, I guess it depends on the pressure you store the hydrogen in the tank;
p_in    = 30                                # Compressor inlet pressure [bar]   # Gabriele: 50 is fine I think, but please check this with Christian. the PEM electrolyzer at Empa works at 30 bar actually
L_is_C  = k/(k-1)*R_H2*T_in_H2*((p_out/p_in)**((k-1)/k)-1) # Specific work compressor [kJ/kg] 

# Maximal TANK energy capacity (J) => 2 days storage capacity
S_TANK_max = E_demand_day * 4 * deltat                  # E_demand_day in [Wh] #Gabriele: E_demand_day is in [Wh], correct? why do you multiply for deltat? 
S_TANK_H2_max = S_TANK_max / HHV                        # equivalent in kg_H2

# Maximal FC size as the maximal power demand divided by eff in [W]
#S_FC_max = P_peak_max / eta["FC"]
S_FC_max = 300*1000  
"""
Range PEMFC = [10W;1MW] from 2021_cigolotti Comprehensive Review on Fuel Cell 
Technology for Stationary Applications
"""
# Prompt the user to choose the scenario
scenario_choice = input("Enter 'grid' for grid-connected scenario or 'off-grid' for an off-grid scenario: ").lower()
# Check the user's choice and set relevant parameters accordingly
if scenario_choice == 'grid':
    P_imp_ub = GRB.INFINITY # Upper bound for P_imp in grid-connected scenario
    P_exp_ub = S_FC_max # Upper bound for P_exp in grid-connected scenario
elif scenario_choice == 'off-grid':
    P_imp_ub = 0  # No imported power in off-grid scenario
    P_exp_ub = 0  # No exported power in off-grid scenario
else:
    print("Invalid scenario choice. Please enter 'grid' or 'off-grid'.")
    sys.exit(1)  # Exit the program in case of an invalid choice
#------------------------------------------------------------------------------
# Define the optimization problem
m = Model()
#------------------------------------------------------------------------------
# Design variables-------------------------------------------------------------
# 1D-Variables
Area_PV = m.addVar(lb=0, ub=Area_PV_max, name='Area_PV')
S_ELY   = m.addVar(lb=0, ub=S_ELY_max,   name='S_ELY')
S_TANK  = m.addVar(lb=0, ub=S_TANK_max,  name='S_TANK')
S_FC    = m.addVar(lb=0, ub=S_FC_max,    name='S_FC')

# Time dependent variables
P_ELY   = m.addVars(nHours, lb=0, ub=S_ELY_max,  name='P_ELY')
E_TANK  = m.addVars(nHours, lb=0, ub=S_TANK_max, name='E_TANK')
P_FC    = m.addVars(nHours, lb=0, ub=S_FC_max,   name='P_FC')    
    
P_imp   = m.addVars(nHours, lb=0, ub=P_imp_ub,   name="P_imp")     # Imported electricity from the Grid in [W]
P_exp   = m.addVars(nHours, lb=0, ub=P_exp_ub,   name="P_exp")     # Exported electricity to the Grid in [W] # Check upper-bound; #gabriele: replace max(P_ELY) with S_ELY. In general, avoid using max() when defining constraints or design variables.
"""
Gabriele:
- In Matlab, we set the ub to the max design value, e.g S_ELY_MAX (or max of 
  demand, which you can extract from the time series), and then we define a 
  constraint like:
  
      m.addConstrs((P_exp[t] <= S_ELY for t in range(nHours)), name="P_exp_upper_bound")
  
- I think we already discussed this, but I forgot :) could you explain me 
  again the rationale to limit the exported power to the ely size? I am in 
  favour of limiting the exported power, I just don't remember why we set 
  that as max value.
"""
#------------------------------------------------------------------------------
# Additional variables definition
# maximum PV size needs to be defined to define the upper bound of the design variable.
P_PV = [irradiance[t] * eta['PV'] * Area_PV for t in range(nHours)]  # Hourly PV power generation [W]
P_PV_peak = 1000 * eta['PV'] * Area_PV  # Peak power for investment cost [W]
"""
Gabriele: do you have a reference for the definition of P_PV_peak you adopted? 
I normally consider the peak power from Standard Test Conditions (STC), for 
which irradiance is fixed at 1000 W/m2. I think it is fairer, otherwise you 
have your peak power depending on sun availability, and given that your unit 
price is CHF/kW, it would not be fair to have an investiment cost depending on 
sun availability (how much units you install instead depends on the sun availability)
"""
S_PV = P_PV_peak

# Hydrogen mass flow
# Nominal mass flow hydrogen [kg/s] #Gabriele: is the flow below kg/s or kg/h?
mdot_H2_nom = (S_ELY * eta["ELY"] * deltat) / HHV                
# Mass flow produced hydrogen [kg/h]                    
mdot_H2 = [(P_ELY[t] * eta['ELY'] * deltat) / HHV for t in range(nHours)] 

# Size (S_C) and operating power (P_C) of the compressor in [kW]
S_C = mdot_H2_nom * L_is_C / eta["C"]  
P_C = [mdot_H2[t] * L_is_C / (eta['C'] * deltat) for t in range(nHours)]

#------------------------------------------------------------------------------
# Constraints
#m.addConstr(E_TANK[0] == S_TANK_max) #Gabriele: this constraint 'conflicts' with line 158, I would use only line 158 and delete line 142

for t in range(1, nHours):                                                     # Gabriele: as you know, vectors in python starts at index 0, so please make sure you are always imposing conditions also for the initial timestep
    # constraint for the energy balance in time over the storage component
    m.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - ((P_FC[t] * deltat) / eta['FC']), name='HESS Balance') #Gabriele: can you please check in the results that this energy balance is respected? Please check over two/three random timesteps, and also over the year, i.e. integral for PEM power, tank energy and FC power. This should differ by the efficiency terms. 
    
    # constraint for the ELY power not to exceed the storage capacity
    m.addConstr(P_ELY[t]  <= (S_TANK - E_TANK[t-1]) / deltat, name='ELY')      # Gabriele: this might be correct as it is already specified in line 148, but please double check that the efficiency is not needed in this inequality
    
    # constraint for the FC to only use hydrogen available in the storage
    m.addConstr(P_FC[t]   <= E_TANK[t-1] / deltat, name='FC')                  # Gabriele: same comment as for line 150

for t in range(nHours):
    m.addConstr((P_ELY[t]  <= S_ELY),  name="upper_Size_Constraint_ELY")
    m.addConstr((E_TANK[t] <= S_TANK), name= "upper_Size_Constraint_TANK")
    m.addConstr((P_FC[t]   <= S_FC ),  name= "upper_Size_Constraint_FC")

# Initializing FC power
m.addConstr(P_FC[0]   <= E_TANK[0] / deltat)                                   # Gabriele: this constraint does not belong to the for loop => removed from for loop

# constraint for H2 storage equal at final and last time step (periodicity)
m.addConstr(E_TANK[0] == E_TANK[nHours-1], name='Periodicity') # 08.02: added name to cosntraint

# Overall energy balance: left => consumers | right => generators
m.addConstrs((P_ELY[t] + P_C[t] + P_exp[t] + P_demand[t] <= P_PV[t] + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance') 

"""
Gabriele: I see here you added a name for the constraints, while you did not 
in other lines. Please make sure that a name is not always needed and that 
constraints are not overwritten if a name is not given.
"""
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Define cost functions
#------------------------------------------------------------------------------

# Installation costs in [€/y]
cost_inst = 0  # Initialize cost_inst outside the loop

# Capital recovery factor
CRF = (annual_interest_rate * (1 + annual_interest_rate)**project_lifetime) / ((1 + annual_interest_rate)**project_lifetime - 1)

for component in ["PV", "ELY", "C", "TANK", "FC"]:
    S_i = globals()["S_" + component]  # Access the variables dynamically
    cost_inst += (S_i * UP[component]) / life[component]   

# Operation costs in [€/y]
cost_op = sum((P_imp[t]/10**6 * cost_imp_el[t] - P_exp[t]/10**6 * cost_exp_el[t]) for t in range(nHours))

# Maintenance cost in [€/y]
cost_maint = 0

# Startup cost in [€/y]
cost_startup = 0

# Total annual costs in [€/y]
cost = cost_inst + cost_op + cost_maint + cost_startup 

#------------------------------------------------------------------------------
# Define the objective function: minimal costs
m.setObjective(cost, GRB.MINIMIZE) 

#------------------------------------------------------------------------------
# Solve optimization problem
m.optimize() 
m.printStats()

"""
Gabriele: overall, the code makes sense to me and I did not find any critical 
mistakes. Please perform these verifications of the results:
   (i) for a simple energy demand and off-grid scenario, compare your estimations 
       (pen and paper, or Excel) and the code results;
  (ii) from the code results in terms of size and operation, calculate from an 
       Excel sheet the total cost and compare it with the objective function 
       calculated in the optimization (they have to be equal)
 (iii) for random timesteps, verify one by one that your constraints are all 
       satisfied. For example, that power balance is correct, that FC power is 
       lower than 'available power' from the tank, etc.
 (iv) verify the energy balances over the whole year. This can be done by 
      integrating the different power flows. Please be careful to consider the 
      efficiencies. 

"""
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
# Retrieve values of variables for further analysis----------------------------

# For Gurobi tupledict object
P_ELY  = [P_ELY[t].X         for t in range(nHours)]
P_FC   = [P_FC[t].X          for t in range(nHours)]
E_TANK = [E_TANK[t].X        for t in range(nHours)]
P_imp  = [P_imp[t].X         for t in range(nHours)]
P_exp  = [P_exp[t].X         for t in range(nHours)]

# For lists with Gurobi LinExpr 
P_PV    = [P_PV[t].getValue()    for t in range(nHours)]
P_C     = [P_C[t].getValue()     for t in range(nHours)]
mdot_H2 = [mdot_H2[t].getValue() for t in range(nHours)]

# For Gurobi var object
Area_PV = Area_PV.X
S_ELY   = S_ELY.X
S_TANK  = S_TANK.X
S_FC    = S_FC.X

# For Gurobi LinExpr
cost        = cost.getValue()
cost_inst   = cost_inst.getValue()
cost_op     = cost_op.getValue()
mdot_H2_nom = mdot_H2_nom.getValue()
P_PV_peak   = P_PV_peak.getValue()
S_C         = S_C.getValue()
S_PV        = S_PV.getValue()
#------------------------------------------------------------------------------
# Calculate the Levelized Cost Of Hydrogen

# LCOH = (cost - cost_startup) / sum(mdot_H2) # From Roxanne 
"""
Definition of LCOH in my model will differ from Roxanne's formulation because
PV can also directly feed the demand or in other words, the demand is not only
satisfied by hydrogen.
"""
#------------------------------------------------------------------------------
# Display results

print("Area_PV = {:.2f} square meters".format(Area_PV))
#------------------------------------------------------------------------------
# Plot main results------------------------------------------------------------

# Import the plotting module
from plotting_module import plot_power_generation, plot_component_sizes, plot_HESS_results, plot_costs_and_prices

# Call the plotting functions as needed
plot_power_generation(P_PV, P_imp, P_exp, df_input)
plot_component_sizes(Area_PV, Area_PV_max, S_ELY, S_FC, S_TANK, S_ELY_max, S_FC_max, S_TANK_max)
plot_HESS_results(P_PV, P_ELY, S_ELY, S_ELY_max, P_FC, S_FC, S_FC_max, E_TANK, S_TANK, S_TANK_max, df_input)
plot_costs_and_prices(cost_inst, cost_op, cost, cost_maint, cost_startup, df_input)

#------------------------------------------------------------------------------
# Export results to excel -----------------------------------------------------

from results_export import export_optimization_results

variable_names = ['irradiance', 'P_demand', 'P_PV', 'P_imp', 'cost_imp_el',
                  'P_exp', 'cost_exp_el', 'P_ELY', 'mdot_H2', 'P_C', 'E_TANK',
                  'P_FC']

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

# Add optimization status to results
results['status'] = m.status

# Define the path to the results directory
results_directory = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\results'

# Export results
export_optimization_results(variable_names, results, results_directory)