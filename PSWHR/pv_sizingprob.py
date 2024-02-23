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
import seaborn as sns # for color palettes



# Introduction: sizing and operational optimization problem
# This model solves a sizing and operational problem with a MILP approach
# The design variable considered is the size of the PV field and the
# hydrogen storage operation.
# The electrolyzer size and battery capacity are instead fixed;

# Set Seaborn style
sns.set(style="whitegrid")
palette = sns.color_palette("husl", 8)  # You can choose a different palette if needed

#------------------------------------------------------------------------------
# Import data
#------------------------------------------------------------------------------
# path to the module 'import_data' => change path if needed
sys.path.append(r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\functions')
from import_data import get_data
df_input, df_demand, irradiance, P_demand = get_data() # Import data and generate variables using the get_data function

#------------------------------------------------------------------------------
# Input parameters
#------------------------------------------------------------------------------
# General parameters
k    = 1.4                               # Ratio cp/cv [-]
R_H2 = 4.1242                            # Individual Gas constant H2 [kJ/kg*K]
HHV  = 39.39*3600*10**3                  # Higher heating value of H2 in [J/kg]

project_lifetime     = 25        # in [y] #Gabriele: do you have a reference for this value by chance?
annual_interest_rate = 0.04      # Discount rate, as encouraged by EU, from Rox


# Electricity prices in [EUR/Wh]
cost_imp_el = df_input['price_Eur_MWh'].values / 10**6 # hourly cost to import 
cost_exp_el = df_input['Price_DayAhed'].values / 10**6 # hourly cost to export 

# Time parameters
nHours = len(df_input)                             # number of hours simulated
Time   = np.arange(1, nHours + 1)                  # time vector
days   = nHours / 24                               # number of days
weeks  = days / 7                                  # number of weeks
deltat = 3600                                      # time step (s)

# Efficiencies of components in [-]
eta = {'PV': 0.21,'ELY': 0+1*0.6, 'C': 0+1*0.7526, 'TANK': 1, 'FC': 0+1*0.5}

# Lifetime of the components in [years]
life = {'PV': 25,'ELY': 20, 'C': 20, 'TANK': 35, 'FC': 20} #Gabriele: do you have references for these values? Can you please check the lifetime values for ely and FC

# Unit prices of components / capital costs in [€/W] => from Roxanne
UP = {
    'PV': 0.8,          # Unit price PV in [€/W]
    'ELY': 1.1,         # Unit price electrolyser [€/W]
    'C': 0.05,          # CHECK FOR DATA
    'TANK': 1*0.1644/HHV, # Unit price hydrogen storage [EUR/J], 1644 [EUR/kgH2]
    'FC': 1,            # CHECK FOR DATA 
    'HP': 0.0576,       # Unit price heat pump
}

# Daily energy demand in [Wh] - P_demand is in [W]
E_demand_day = sum(P_demand)/(4*days) #Gabriele: this is [Wh], correct?
# Maximal daily power demand in [W]
P_peak_max = np.max(P_demand)

# Defining maximal SIZES of the components-------------------------------------
# maximum PV size needs to be defined to define the upper bound of the design variable.
# Area_PV_max  = P_peak_max / (eta['PV'] * np.mean(irradiance)) # Maximum PV area [m2] 
Area_PV_max  = 0.005*100000 # Maximum PV area [m2]

# Electrolyser maximal nominal power (W) very large (not realistic)
S_ELY_max = 500*1000                         # Maximal Electrolyzer size in [W]

# Calculating the spezific work of the compressor, from Minutillo et al. 2021
# (Analyzing the levelized cost of hydrogen eq 1+2) => from Roxanne 
T_in_H2 = 70 + 273.15                       # H2 temperature (=T_cat=T_an) [K] #Gabriele: this could probably be 65 degC, see for example: https://www.sciencedirect.com/science/article/pii/S0360319923015410
p_out   = 820                               # Compressor outlet pressure [bar] #Gabriele: is your p_out at 820 bar? Please see in the litearture what value is set here, I guess it depends on the pressure you store the hydrogen in the tank;
p_in    = 50                                # Compressor inlet pressure [bar] #Gabriele: 50 is fine I think, but please check this with Christian. the PEM electrolyzer at Empa works at 30 bar actually
L_is_C  = k/(k-1)*R_H2*T_in_H2*((p_out/p_in)**((k-1)/k)-1) # Specific work compressor [kJ/kg] 

# Maximal TANK energy capacity (J) => 2 days storage capacity
S_TANK_max = 10000*E_demand_day * 4 * deltat                  # E_demand_day in [Wh] #Gabriele: E_demand_day is in [Wh], correct? why do you multiply for deltat?
S_TANK_H2_max = S_TANK_max / HHV                        # equivalent in kg_H2

# Maximal FC size as the maximal power demand divided by eff in [W]
#S_FC_max = P_peak_max / eta["FC"]
S_FC_max = 300*1000

#------------------------------------------------------------------------------
# Define the optimization problem
PV_sizingprob = Model()
#------------------------------------------------------------------------------
# Design variables-------------------------------------------------------------
# 1D-Variables
Area_PV = PV_sizingprob.addVar(lb=0, ub=Area_PV_max, name='Area_PV')
S_ELY   = PV_sizingprob.addVar(lb=0, ub=S_ELY_max,   name='S_ELY')
S_TANK  = PV_sizingprob.addVar(lb=0, ub=S_TANK_max,  name='S_TANK')
S_FC    = PV_sizingprob.addVar(lb=0, ub=S_FC_max,    name='S_FC')

# Time dependent variables
P_ELY   = PV_sizingprob.addVars(nHours, lb=0, ub=S_ELY_max,    name='P_ELY')
E_TANK  = PV_sizingprob.addVars(nHours, lb=0, ub=S_TANK_max,   name='E_TANK')
P_FC    = PV_sizingprob.addVars(nHours, lb=0, ub=S_FC_max,     name='P_FC')    
P_imp   = PV_sizingprob.addVars(nHours, lb=0, ub=1*1000*1000 + 0*GRB.INFINITY, name="P_imp")   # Imported electricity from the Grid in [W]
P_exp   = PV_sizingprob.addVars(nHours, lb=0, ub=1*1000*1000 + 0*max(P_ELY),   name="P_exp")   # Check upper-bound; #gabriele: replace max(P_ELY) with S_ELY. In general, avoid using max() when defining constraints or design variables.
#Gabriele: In Matlab, we set the ub to the max design value, e.g S_ELY_MAX (or max of demand, which you can extract from the time series), and then we define a constraint like:
#Gabriele: PV_sizingprob.addConstrs((P_exp[t] <= S_ELY for t in range(nHours)), name="P_exp_upper_bound")
#Gabriele: I think we already discussed this, but I forgot :) could you explain me again the rationale to limit the exported power to the ely size? I am in favour of limiting the exported power, I just don't remember why we set that as max value

#------------------------------------------------------------------------------
# Additional variables definition
# maximum PV size needs to be defined to define the upper bound of the design variable.
P_PV = [irradiance[t] * eta['PV'] * Area_PV for t in range(nHours)]  # Hourly PV power generation [W]
P_PV_peak = max(irradiance) * eta['PV'] * Area_PV  # Peak power for investment cost [W] #Gabriele: in optimization problems, try to avoid using the max() function. Besides, the peak power 
#Gabriele: do you have a reference for the definition of P_PV_peak you adopted? I normally consider the peak power from Standard Test Conditions (STC), for which irradiance is fixed at 1000 W/m2. I think it is fairer, otherwise you have your peak power depending on sun availability, and given that your unit price is CHF/kW, it would not be fair to have an investiment cost depending on sun availability (how much units you install instead depends on the sun availability)
S_PV = P_PV_peak


# Hydrogen mass flow
# Nominal mass flow hydrogen [kg/s] #Gabriele: is the flow below kg/s or kg/h?
mdot_H2_nom = (S_ELY * eta["ELY"] * deltat) / HHV     # wird nicht ausgerechnet       
# Mass flow produced hydrogen [kg/h]                    
mdot_H2 = [(P_ELY[t] * eta['ELY'] * deltat) / HHV for t in range(nHours)] 

# Size (S_C) and operating power (P_C) of the compressor in [kW]
S_C = mdot_H2_nom * L_is_C / eta["C"]  
P_C = [mdot_H2[t] * L_is_C / (eta['C'] * deltat) for t in range(nHours)]

#------------------------------------------------------------------------------
# Constraints
PV_sizingprob.addConstr(E_TANK[0] == 0) #Gabriele: this constraint 'conflicts' with line 158, I would use only line 158 and delete line 142

for t in range(1, nHours): #Gabriele: as you know, vectors in python starts at index 0, so please make sure you are always imposing conditions also for the initial timestep
    # constraint for the energy balance in time over the storage component
    PV_sizingprob.addConstr(E_TANK[t] == E_TANK[t-1] + P_ELY[t] * eta['ELY'] * deltat - P_FC[t] * deltat / eta['FC']) #Gabriele: can you please check in the results that this energy balance is respected? Please check over two/three random timesteps, and also over the year, i.e. integral for PEM power, tank energy and FC power. This should differ by the efficiency terms. 
    # constraint for the ELY power not to exceed the storage capacity
    PV_sizingprob.addConstr(P_ELY[t]  <= (S_TANK - E_TANK[t-1]) / deltat) #Gabriele: this might be correct as it is already specified in line 148, but please double check that the efficiency is not needed in this inequality
    # constraint for the FC to only use hydrogen available in the storage
    PV_sizingprob.addConstr(P_FC[t]   <= E_TANK[t-1] / deltat) #Gabriele: same comment as for line 150
    PV_sizingprob.addConstr(P_FC[0]   <= E_TANK[0] / deltat) #Gabriele: this constraint does not belong to the for loop
    
    PV_sizingprob.addConstr((P_ELY[t]  <= S_ELY),  "upper_Size_Constraint_ELY")
    PV_sizingprob.addConstr((E_TANK[t] <= S_TANK), "upper_Size_Constraint_TANK")
    PV_sizingprob.addConstr((P_FC[t]   <= S_FC ),   "upper_Size_Constraint_FC")
    #PV_sizingprob.addConstr((P_imp[t]   <= 10000.0 ),   "upper_Size_Constraint_Grid")

# constraint for H2 storage equal at final and last time step (periodicity)
PV_sizingprob.addConstr(E_TANK[0] == E_TANK[nHours-1]) 

# Overall energy balance: left => consumers | right => generators
PV_sizingprob.addConstrs((P_ELY[t] + P_C[t] + P_exp[t] + P_demand[t] <= P_PV[t] + P_imp[t] + P_FC[t]  for t in range(nHours)), name='EnergyBalance') #Gabriele: I see here you added a name for the constraints, while you did not in other lines. Please make sure that a name is not always needed and that constraints are not overwritten if a name is not given.

#------------------------------------------------------------------------------
# Define cost functions
#------------------------------------------------------------------------------
components = ["PV", "ELY", "C", "TANK", "FC"]

# Installation costs in [€/y]
cost_inst = 0  # Initialize cost_inst outside the loop

for component in components:
    S_i = globals()["S_" + component]  # Access the variables dynamically
    cost_inst += (S_i * UP[component]) / life[component]   

cost_inst = 0

# Operation costs in [€/y]
cost_op = 0
#cost_op = sum((P_imp[t] * cost_imp_el[t] - P_exp[t] * cost_exp_el[t]) for t in range(nHours))
#cost_op = sum((P_imp[t] * cost_imp_el[t] - 0*P_exp[t] * cost_exp_el[t]) for t in range(nHours))
#cost_op = cost_op + 12*max(P_imp)*16.0 # max() bad style by PC54!
'''
P_max_imp = 0
for t in range(nHours):
    # if P_imp[t] > P_max_imp:
        # P_max_imp = P_imp[t]
        PV_sizingprob.addConstr((P_max_imp <= P_imp[t]) >> (P_max_imp==P_imp[t]), "max_grid_power_constraint")

cost_op = cost_op + 12*P_max_imp*16.0
'''


# Maintenance cost 
# cost_maint = 

# Total annual costs in [k€/y]
cost = cost_inst + cost_op #+ cost_maint + cost_startup

#------------------------------------------------------------------------------
# Define the objective function: minimal costs
PV_sizingprob.setObjective(cost, GRB.MINIMIZE) 

#------------------------------------------------------------------------------
# Solve optimization problem
PV_sizingprob.optimize() 

#Gabriele: overall, the code makes sense to me and I did not find any critical mistakes. Please perform these verifications of the results:
#Gabriele: (i) for a simple energy demand and off-grid scenario, compare your estimations (pen and paper, or Excel) and the code results;
#Gabriele: (ii) from the code results in terms of size and operation, calculate from an Excel sheet the total cost and compare it with the objective function calculated in the optimization (they have to be equal)
#Gabriele: (iii) for random timesteps, verify one by one that your constraints are all satisfied. For example, that power balance is correct, that FC power is lower than 'available power' from the tank, etc.
#Gabriele: (iv) verify the energy balances over the whole year. This can be done by integrating the different power flows. Please be careful to consider the efficiencies. 


#------------------------------------------------------------------------------
# Check the optimization status
if PV_sizingprob.status == GRB.OPTIMAL:
    print("Optimal solution found.")
elif PV_sizingprob.status == GRB.INFEASIBLE:
    print("The model is infeasible.")
elif PV_sizingprob.status == GRB.UNBOUNDED:
    print("The model is unbounded.")
elif PV_sizingprob.status == GRB.TIME_LIMIT:
    print("Time limit reached.")

# Print the status message for more details => check here: https://www.gurobi.com/documentation/current/refman/optimization_status_codes.html 
print(f"Optimization status: {PV_sizingprob.status}")
#------------------------------------------------------------------------------
#Display results

all_vars = PV_sizingprob.getVars()
values   = PV_sizingprob.getAttr("X", all_vars)
names    = PV_sizingprob.getAttr("VarName", all_vars)

results = {}  # Create a dictionary to store the results

for name, val in zip(names, values):
    results[name] = val
    #print(f"{name} = {val}")

print("Area_PV = {:.2f} square meters".format(Area_PV.X))
print("Max grid power peak = {:.2f} kW".format(max(P_imp)))

try:
    #------------------------------------------------------------------------------
    # Plot main results------------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(20, 10))

    # Plot 1: Size and operation of the ELY
    ax1.plot(range(len(df_input)), [P_ELY[t].X / 1000 for t in range(len(df_input))], label='ELY Electrolyzer (kW)', color=palette[1])
    ax1.plot(range(len(df_input)), [S_ELY.X / 1000] * len(df_input), label='ELY size (kW)', linestyle='--', color=palette[3])
    ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Power [kW]')
    ax1.set_ylim([0, S_ELY.X * 1.1 / 1000])
    ax1.legend(loc='upper right')

    # Plot 1: Size and operation of the FC
    ax2.plot(range(len(df_input)), [P_FC[t].X / 1000 for t in range(len(df_input))], label='Fuel Cell (kW)', color=palette[0])
    ax2.plot(range(len(df_input)), [S_FC.X / 1000] * len(df_input), label='FC size (kW)', linestyle=':', color=palette[2])
    ax2.set_xlabel('Time [h]')
    ax2.set_ylabel('Power [kW]')
    ax2.set_ylim([0, S_FC.X * 1.1 / 1000])
    ax2.legend(loc='upper right')

    # Plot 3: Energy and TANK Size
    ax3.plot(range(len(df_input)), [E_TANK[t].X / 3600 / 1000 for t in range(len(df_input))], label='TANK (kWh)', color=palette[5])
    ax3.plot(range(len(df_input)), [S_TANK.X / 1000 / 3600] * len(df_input), label='TANK size (kWh)', linestyle='-.', color=palette[6])
    ax3.set_xlabel('Time [h]')
    ax3.set_ylabel('Energy [kWh]')
    ax3.set_ylim([0, S_TANK.X * 1.1 / 1000 / 3600])
    ax3.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    # Sizes of the components------------------------------------------------------
    # Assuming 'S_ELY', 'S_TANK', 'S_FC' are Gurobi variables
    S_ELY_value = S_ELY.X
    S_TANK_value = S_TANK.X
    S_FC_value = S_FC.X

    # Create a bar plot
    components = ['S_ELY', 'S_TANK', 'S_FC']
    component_sizes = [S_ELY_value, S_TANK_value, S_FC_value]

    # Bar plot
    plt.figure(figsize=(6, 4))
    plt.bar(components, component_sizes, color=['blue', 'green', 'orange', 'red'])
    plt.title('Sizes of Components')
    plt.xlabel('Components')
    plt.ylabel('Size')
    plt.show()

    import matplotlib.pyplot as plt

    # Plot visualizing the electricity generated, imported, and exported-----------
    # Retrieve the values of P_PV, P_imp, and P_exp => stacked plot here
    P_PV_values = [P_PV[t].getValue()/1000 for t in range(nHours)]
    P_imp_values = [P_imp[t].X/1000 for t in range(nHours)]
    P_exp_values = [P_exp[t].X/1000 for t in range(nHours)]

    # Create subplots
    fig, axs = plt.subplots(2, 1, figsize=(20, 10), sharex=True)

    # Plot P_PV and P_imp in the first subplot
    axs[0].plot(range(len(df_input)), P_PV_values, label='PV Generation (W)', color='orange')
    axs[0].plot(range(len(df_input)), P_imp_values, label='Imported Power (W)', color='blue')
    axs[0].set_ylabel('Power [kW]')
    axs[0].legend(loc='upper left')
    axs[0].set_title('PV Generation and Imported Power')

    # Plot P_exp in the second subplot
    axs[1].plot(range(len(df_input)), P_exp_values, label='Exported Power (W)', color=palette[5])
    axs[1].set_xlabel('Time [h]')
    axs[1].set_ylabel('Power [kW]')
    axs[1].legend(loc='upper left')
    axs[1].set_title('Exported Power')

    plt.show()


    # Plot costs and electricity prices--------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10))
    # Plot 1: Import and Export Prices
    ax1.plot(range(len(df_input)), df_input['price_Eur_MWh'], label='Import Price (EUR/MWh)', color=palette[3])
    ax1.plot(range(len(df_input)), df_input['Price_DayAhed'], label='Export Price (EUR/MWh)', color=palette[4])
    ax1.set_xlabel('Time [h]')
    ax1.set_ylabel('Price [EUR/Wh]')
    ax1.legend()

    # Plot 2: Costs Overview (Bar Plot)
    bar_width = 0.4

    bar1 = ax2.bar(0, cost_inst.getValue()/1000, bar_width, label='Installation Cost (k€/y)', color=palette[0])
    bar2 = ax2.bar(1, cost_op.getValue()/1000, bar_width, label='Operation Cost (k€/y)', color=palette[1])
    bar3 = ax2.bar(2, cost.getValue()/1000, bar_width, label='Total Cost (k€/y)', color=palette[2])

    ax2.set_ylabel('Cost [k€/y]')
    ax2.legend()
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(['Installation Cost', 'Operation Cost', 'Total Cost'])
    plt.show()
except:
    print("An error occured during plotting results")
#------------------------------------------------------------------------------
# Export results to excel -----------------------------------------------------

from results_export import export_optimization_results

# Extract the values of decision variables and other variables
variable_names = ['irradiance', 'P_demand', 'P_PV', 'P_imp', 'cost_imp_el', 
                  'P_exp', 'cost_exp_el','P_ELY', 'mdot_H2', 'P_C', 'E_TANK', 
                  'P_FC', 'status']

variable_values = {name: [] for name in variable_names}

for t in range(nHours):
    variable_values['irradiance'].append(irradiance[t])
    variable_values['P_demand'].append(P_demand[t])
    variable_values['P_PV'].append(P_PV[t].getValue())
    variable_values['P_imp'].append(P_imp[t].X)
    variable_values['cost_imp_el'].append(cost_imp_el[t])
    variable_values['P_exp'].append(P_exp[t].X)
    variable_values['cost_exp_el'].append(cost_exp_el[t])
    variable_values['P_ELY'].append(P_ELY[t].X)
    variable_values['mdot_H2'].append(mdot_H2[t])
    variable_values['P_C'].append(P_C[t])
    variable_values['E_TANK'].append(E_TANK[t].X)
    variable_values['P_FC'].append(P_FC[t].X)

# Add optimization status to variable_values
variable_values['status'] = PV_sizingprob.status

# Define the path to the results directory
results_directory = r'C:\Users\peter_c\Desktop\test\results'

# Export results
export_optimization_results(variable_names, variable_values, results_directory)
