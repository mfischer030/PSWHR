import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import numpy as np

'''
Comments Gabriele:
The code is written in kW and kWh
The code works well if short period of times are considered, e.g. 1 week
For 1-year analyis, the solution convergence if battery size is fixed, i.e. battery dispatch problem. 
I am not sure the code coverges for 1-year analysis and sizing+operation optimization, the simulation was taking too long. This might be to the adopted random demand, to the parameters (like the random costs I used), to the way the code is written (can be improved?), and/or to the big-M constraint, which can be a source of instability (https://www.gurobi.com/documentation/current/refman/dealing_with_big_m_constra.html)
Nonetheless, I suspect that if this method is embedded in a larger model, e.g. Maxime's case of the hydrogen storage system, it might need larger computational power than a regular laptop.

Note: the battery efficiencies might be counted twice at some point in the code, please have a look
'''

plt.close('all')

# Parameters
max_sim_time  = 2*60  # the 1-year test can take longer than this time for the boundary conditions and parameter set in the initial version of the code
hours_in_year = 8760
eff_b_ch      = 0.95  # Battery charging efficiency
eff_b_disch   = 0.95  # Battery discharging efficiency
eff_b_sd      = 0.99  # Battery self-discharge efficiency
c_imp_low     = 0.1   # Cost in CHF/kWh for low grid usage - very low to test model
c_imp_high    = 5.0   # Cost in CHF/kWh for high grid usage - very high to test model
epsilon       = 1e-4  # Small value to avoid numerical issues
P_dem_peak    = 3     # here we set the peak for the randonm demand (not needed if demand is imported)
C_b_max       = 50    # Max battery capacity in kWh y
C_b_min       = 0     # Min battery capacity in kWh - typically set to 0, but if set equal to C_b_max can be used to impose a fixed battery capacity without modifying the code 
UP_b          = 500   # unit price of the battery in CHF/kWh - I selected a low price to promote battery installation in the optimization
lifetime_b    = 10    # lifetime of the battery in years - I selected a large lifetime to promote battery installation
SOC_max       = 0.8   # max state of charge of the battery to increase lifetime
SOC_min       = 0.2   # min state of charge of the battery to increase lifetime
threshold_hours = 3500  # Threshold for low/high grid usage in hours

# these two coefficients are needed for the big-M constraints
# Tipically, the M value is selected as a large value, but too large values must be avoided to avoid numerical instabilities
# the selection of the M values depend on the problem you are tackling
M_hours     = 2*P_dem_peak    # used for the big-M constraint on the grid_usage
M_threshold = hours_in_year   # used for the big-M for the threshold

# reduced test - I reduce vector length to visualize better results and have faster tests
# factor = 1                            # set factor to 1 to run the 1-year test
# factor = 8760//(24*30*3)              # factor to reduce analysis to 1 quarter
factor = 8760//(24*30)                  # factor to reduce analysis to 1 month
# factor = 8760//(24*7)                 # factor to reduce analysis to 1 week
# factor = 8760//(24)                   # factor to reduce analysis to 1 day

hours_in_year   = hours_in_year // factor
threshold_hours = threshold_hours // factor
lifetime_b=lifetime_b * factor

# Random demand
np.random.seed(0)
E_dem = np.random.uniform(low=0, high=P_dem_peak, size=hours_in_year)  # Random demand in kWh

# Optimization model
model = gp.Model("BatteryDispatch")
model.setParam('TimeLimit', max_sim_time)

# Decision variables-----------------------------------------------------------

P_imp   = model.addVars(hours_in_year, name="P_imp")                           # Power imported from the grid
P_ch    = model.addVars(hours_in_year, name="P_ch")                            # Power charged to the battery
P_disch = model.addVars(hours_in_year, name="P_disch")                         # Power discharged from the battery
E_b     = model.addVars(hours_in_year + 1, lb=0, ub=C_b_max, name="E_b")       # Energy in the battery at each time step
C_b     = model.addVar(lb=0, ub=C_b_max, name="C_b")                           # Battery Capacity

grid_usage = model.addVars(hours_in_year, vtype=GRB.BINARY, name="grid_usage") # Binary variable for grid usage indicator, 1 if used, i.e. P_imp>0, 0 otherwise
h_usage    = model.addVar(lb=0, ub=hours_in_year, name="h_usage")              # Total hours of grid usage
high_usage = model.addVar(vtype=GRB.BINARY, name="high_usage")                 # Binary variable indicating whether grid usage is above the threshold, 1 if above, 0 if below

# Constraints------------------------------------------------------------------

# First we model the physical system

# Energy balance for the battery storage
model.addConstrs(E_b[i + 1] == E_b[i]* eff_b_sd + P_ch[i] * eff_b_ch - P_disch[i] / eff_b_disch for i in range(hours_in_year)) 

# Constraint for discharged power lower than available power
model.addConstrs(P_disch[i] <= E_b[i] for i in range(hours_in_year)) 

# The following constraints increase computational time significantly I think, I guess in particular the SOC_min one
model.addConstrs(E_b[i] <= SOC_max * C_b for i in range(hours_in_year)) 
model.addConstrs(E_b[i] >= SOC_min * C_b for i in range(hours_in_year))

# Energy balance
model.addConstrs(P_disch[i] / eff_b_disch + P_imp[i] >= E_dem[i] + P_ch[i]  for i in range(hours_in_year)) 

# Periodicity 
model.addConstr(E_b[hours_in_year] == E_b[0]) 




# here we add the constraints needed for the grid usage tarif selection
model.addConstr(h_usage == gp.quicksum(grid_usage[i] for i in range(hours_in_year))) # not sure this has to be a constraint, could be a normal calculation maybe 

# we use the big-M approach to define the binary variable grid_usage
for i in range(hours_in_year):
    model.addConstr(P_imp[i] <= M_hours * grid_usage[i])                 # here we constraint grid_usage to take value 1 if P_imp is larger than zero. If P_imp is positive, this constraint is satisfied only if grid_usage is 1
    model.addConstr(P_imp[i] >= epsilon - M_hours * (1 - grid_usage[i])) # here we constraint grid_usage to 0 if P_imp is zero. 
# Here I choose M_hours in the order of the peak demand

# defining the value of the binary variable high-usage based on the hours of operation, h_usage
model.addConstr(h_usage - threshold_hours <= M_threshold * high_usage)# if h_usage is larger than threshold, the inequality is respected only if high_usage takes value 1. 
# here I choose the M_threshold in the order of the hours in the year

# Objective function-----------------------------------------------------------
cost_inv = C_b*UP_b/lifetime_b
cost_op  = gp.quicksum(P_imp[i] * (c_imp_low * (1 - high_usage) + c_imp_high * high_usage) for i in range(hours_in_year)) # here we define the import cost as a function of the high_usage binary variable
cost     = cost_inv+cost_op

model.setObjective(cost, GRB.MINIMIZE)

# Solve the model
model.optimize()
    
if model.status == gp.GRB.OPTIMAL:
    print("Optimal solution found.")
    print(f"Total cost: {model.objVal} CHF")
    print(f"Total grid usage hours: {h_usage.X}")
elif model.status == gp.GRB.TIME_LIMIT:
    print("Time limit reached. Best solution found:")
    print(f"Total cost: {model.objVal} CHF")
    print(f"Total grid usage hours: {h_usage.X}")
else:
    print("Solver ended with status:", model.status)

#######################################################
## Plotting ##

# Post-optimization: Extracting values for plotting
energy_demand = [E_dem[i] for i in range(hours_in_year)]
power_supplied = [P_disch[i].X + P_imp[i].X - P_ch[i].X for i in range(hours_in_year)]
import_power = [P_imp[i].X for i in range(hours_in_year)]
charging_power = [P_ch[i].X for i in range(hours_in_year)]
discharging_power = [P_disch[i].X for i in range(hours_in_year)]
battery_energy = [E_b[i].X for i in range(hours_in_year + 1)]
selected_electricity_cost = [c_imp_low if h_usage.X <= 3500 else c_imp_high for _ in range(hours_in_year)]
grid_usage_binary = [grid_usage[i].X for i in range(hours_in_year)]
cumulative_hours_usage = np.cumsum(grid_usage_binary)
selected_battery_capacity = C_b.X  # Get the optimized battery capacity
cost_inv = C_b.X * UP_b / lifetime_b
cost_op = model.getObjective().getValue() - cost_inv

# First set of tile plots
fig, axs = plt.subplots(3, 1, figsize=(12, 18))

# Tile 1: Energy demand and sum of discharged battery power and import grid vs time
axs[0].plot(energy_demand, label='Energy Demand (kWh)', color='blue')
axs[0].plot(power_supplied, label='Sum of Discharged and Imported Power minus charged (kWh)', color='red')
#axs[0].set_title('Energy Demand and Supply Over Time - They should be equal; I guess I counted battery eff twice somewhere')
axs[0].set_xlabel('Time (hours)')
axs[0].set_ylabel('Energy (kWh)')
axs[0].legend()

# Tile 2: Import power, charging power, and discharging power vs time
axs[1].plot(import_power, label='Import Power (kW)', color='green')
axs[1].plot(charging_power, label='Charging Power (kW)', color='orange')
axs[1].plot(discharging_power, label='Discharging Power (kW)', color='purple')
#axs[1].set_title('Power Flows Over Time')
axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('Power (kW)')
axs[1].legend()

# Tile 3: Energy in the battery vs time
axs[2].plot(battery_energy, label='Battery Energy (kWh)', color='cyan')
#axs[2].set_title('Battery State of Charge Over Time')
axs[2].axhline(y=SOC_min * selected_battery_capacity, color='orange', linestyle='--', 
               label=f'Min SOC ({SOC_min * 100}%) Capacity: {SOC_min * selected_battery_capacity:.2f} kWh')
axs[2].axhline(y=SOC_max * selected_battery_capacity, color='green', linestyle='--', 
               label=f'Max SOC ({SOC_max * 100}%) Capacity: {SOC_max * selected_battery_capacity:.2f} kWh')
axs[2].axhline(y=selected_battery_capacity, color='red', linestyle='-', 
               label=f'Selected Capacity: {selected_battery_capacity:.2f} kWh')
axs[2].set_xlabel('Time (hours)')
axs[2].set_ylabel('Energy (kWh)')
axs[2].set_ylim([0, selected_battery_capacity * 1.1])
axs[2].legend()

plt.tight_layout()
plt.show()

# Second set of tile plots
fig, axs = plt.subplots(4, 1, figsize=(12, 24))

# Tile 1: Selected electricity cost vs time
axs[0].plot(selected_electricity_cost, label='Selected Electricity Cost (CHF/kWh)', color='black')
#axs[0].set_title('Electricity Cost Over Time')
axs[0].set_xlabel('Time (hours)')
axs[0].set_ylabel('Cost (CHF/kWh)')
axs[0].legend()

# Tile 2: Electricity imported vs time
axs[1].plot(import_power, label='Imported Electricity (kW)', color='green')
#axs[1].set_title('Imported Electricity Over Time')
axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('Electricity (kW)')
axs[1].legend()

# Tile 3: Grid usage vs time
axs[2].step(range(hours_in_year), grid_usage_binary, label='Grid Usage (Binary)', where='post', color='red')
#axs[2].set_title('Grid Usage Status Over Time')
axs[2].set_xlabel('Time (hours)')
axs[2].set_ylabel('Grid Usage (Binary)')
axs[2].legend()
axs[2].set_ylim([0, 2])


# Tile 4: Cumulative hours of usage vs time
axs[3].plot(cumulative_hours_usage, label='Cumulative Grid Usage Hours', color='blue')
axs[3].axhline(y=threshold_hours, color='red', linestyle='--', label=f'Threshold Hours ({threshold_hours})')
#axs[3].set_title('Cumulative Grid Usage Hours Over Time')
axs[3].set_xlabel('Time (hours)')
axs[3].set_ylabel('Cumulative Hours')
axs[3].set_ylim([0, hours_in_year + (hours_in_year * 0.1)])  # Adjust the y-axis limit to ensure visibility of all lines
axs[3].legend()

plt.tight_layout()
plt.show()

# bar plots for cost
cost_categories = ['Investment Cost', 'Operational Cost']
cost_values = [cost_inv, cost_op]

plt.figure(figsize=(8, 6))
bar_positions = range(len(cost_categories))
plt.bar(bar_positions, cost_values, color=['skyblue', 'lightgreen'])
plt.title('Investment vs Operational Costs')
plt.ylabel('Cost in CHF')
plt.xticks(bar_positions, cost_categories)

for i, value in enumerate(cost_values):
    plt.text(i, value + 0.05 * max(cost_values), f'{value:.2f}', ha='center')

plt.tight_layout()
plt.show()