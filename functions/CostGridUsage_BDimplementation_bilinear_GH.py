import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import numpy as np

# Close all previous plots
plt.close('all')

# Parameters
max_sim_time   = 1*60
months         = 12
hours_in_month = 24*30                  # I define months all of the same length to simplify this model
hours_in_year  = months*hours_in_month  # hours in a year if 12 months of 30 days
threshold_BD   = 3500                   # Threshold for low/high grid usage (BD threshold)
epsilon        = 1e-4
P_dem_peak     = 3                      # with Peak 300 kW, En demand results as 1284245 kWh, which is BD 4280 for peak at 300 
P_imp_ub       = P_dem_peak*24
lifetime_b     = 10
SOC_max        = 1.0
SOC_min        = 0.0
M_hours        = 5 * P_dem_peak
M_threshold    = hours_in_year
M_BD           = 10000                 # Large constant, value could be reduced???
penalty_factor = 10.0                  # penalty for the P_max

# Baseline parameters
"""
eff_b_ch = 0.95
eff_b_disch = 0.95
eff_b_sd = 0.99
c_imp_low = 0.1 # price when BD < thereshold
c_imp_high = 5.0 # price when BD > thereshold
C_b_max = P_dem_peak*48
UP_b = 500
C_b_min = 0
"""

# Adotped parameters
eff_b_ch    = 1.0
eff_b_disch = 1.0
eff_b_sd    = 1.0
c_imp_low   = 0.25  # import price when BD < thereshold
c_imp_high  = 0.2   # import price when BD > thereshold [CHF/kWh]
c_peak_high = 15    # price for peak power when BD > thereshold [CHF/kW/month]
c_peak_low  = 8     # price for peak power when BD > thereshold [CHF/kW/month]
C_b_max     = P_dem_peak*48*5
UP_b        = 500
C_b_min     = 0


# Test scale factor
factor          = 1 # keep it to 1, model not adapted to different factors
hours_in_year //= factor
threshold_BD  //= factor
lifetime_b     *= factor

# Random demand setup
np.random.seed(0)
E_dem = np.random.uniform(low=0, high=P_dem_peak, size=hours_in_year)

# Create the optimization model
model = gp.Model("BatteryDispatch")
model.setParam('TimeLimit', max_sim_time)

# Decision variables
P_imp      = model.addVars(hours_in_year,                       name="P_imp")
P_max_imp  = model.addVars(months, lb=0, ub=P_imp_ub,           name="P_max_imp")
P_ch       = model.addVars(hours_in_year,                       name="P_ch")
P_disch    = model.addVars(hours_in_year,                       name="P_disch")
E_b        = model.addVars(hours_in_year + 1, lb=0, ub=C_b_max, name="E_b")
C_b        = model.addVar(lb=C_b_min, ub=C_b_max,               name="C_b")
high_usage = model.addVar(vtype=GRB.BINARY,                     name="high_usage")

# Constraints for battery dynamics
model.addConstrs(E_b[i + 1] == E_b[i] * eff_b_sd + P_ch[i] * eff_b_ch - P_disch[i] / eff_b_disch for i in range(hours_in_year))
model.addConstrs(P_disch[i] <= E_b[i] for i in range(hours_in_year))
# model.addConstrs(P_ch[i] <= E_b[i] for i in range(hours_in_year))
model.addConstrs(E_b[i] <= SOC_max * C_b for i in range(hours_in_year))
model.addConstrs(E_b[i] >= SOC_min * C_b for i in range(hours_in_year))
# Does inequality constraint lead to artificial P_imp to reduce BD? It seems no
model.addConstrs(P_disch[i] / eff_b_disch + P_imp[i] >= E_dem[i] + P_ch[i] for i in range(hours_in_year))
#model.addConstrs(P_disch[i] / eff_b_disch + P_imp[i] == E_dem[i] + P_ch[i] for i in range(hours_in_year))

model.addConstr(E_b[hours_in_year] == E_b[0])

# Monthly maximum constraint - @ Maxime, this is something you already implemented in your code. In your code it is more accurate as here I assume 12 months of 30 days each
# I was having some issues here with the P_max_imp value selected by the optimizer
# basically, in order to reduce the BD, the optimizer can choose a larger value than the actual one
# for example, for a demand of 10 kWh, I might have a average peak of power of 1kW and so a BD of 10, but the optimizer could choose a larger P_max_imp, e.g. 2 kW, to reduce BD to 5.
# I tried to think of some ways to avoid this but I was not succesfull...
# However, I think the issue is not there if the cost for the P_peak is sufficiently high.
# Anyway, It is important to always verify the results in the post-processing!

for month in range(months):
    start = month * hours_in_month
    end = start + hours_in_month
    for t in range(start, end):
        model.addConstr(P_imp[t] <= P_max_imp[month], name=f"max_monthly_constraint_{month}_{t}")

# imported power and average of imported peak
total_imported_power = gp.quicksum(P_imp[i] for i in range(hours_in_year))
avg_monthly_max = gp.quicksum(P_max_imp[j] for j in range(months)) / months

# BD constraints using binary variable and bigM
model.addConstr(total_imported_power - threshold_BD * avg_monthly_max >= -M_BD * (1 - high_usage), "activate_high_usage")
model.addConstr(total_imported_power - threshold_BD * avg_monthly_max <= M_BD * high_usage, "deactivate_high_usage")

# Objective function with new cost logic

# Modify the objective function

cost_inv = C_b * UP_b / lifetime_b

# This is not linear 
cost_op = gp.quicksum(P_imp[i] * (c_imp_low * (1 - high_usage) + c_imp_high * high_usage) for i in range(hours_in_year))

# This is not linear 
cost_peak= gp.quicksum(P_max_imp[month] * (c_peak_low * (1 - high_usage) + c_peak_high * high_usage) for month in range(months))
model.setObjective(cost_inv + cost_op + cost_peak, GRB.MINIMIZE)

# Solve the model
model.optimize()

"""
# Check results
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
"""

#%% Post-optimization: Extracting values for plotting
monthly_max_demand = []

# Loop through each month
for month in range(months):
    start_index = month * hours_in_month
    end_index = start_index + hours_in_month
    # Make sure to limit the end_index to the length of E_dem to avoid index out of range error
    month_demands = E_dem[start_index:min(end_index, len(E_dem))]
    # Append the maximum demand of this month to the list
    monthly_max_demand.append(max(month_demands))

average_max_demand = sum(monthly_max_demand) / len(monthly_max_demand)
energy_demand = [E_dem[i] for i in range(hours_in_year)]
energy_demand_year= sum(energy_demand)

import_max_opt = [P_max_imp[month].X for month in range(months)]
avg_import_max_opt = np.mean(import_max_opt)
BD_opt=energy_demand_year/avg_import_max_opt

power_supplied = [P_disch[i].X + P_imp[i].X - P_ch[i].X for i in range(hours_in_year)]
import_power = [P_imp[i].X for i in range(hours_in_year)]
energy_imported_year=sum(import_power)

import_max_postprocessing= []

for month in range(months):
    start_index = month * hours_in_month
    end_index = start_index + hours_in_month
    # Make sure to limit the end_index to the length of E_dem to avoid index out of range error
    import_postprocessing = import_power[start_index:min(end_index, len(import_power))]
    # Append the maximum demand of this month to the list
    import_max_postprocessing.append(max(import_postprocessing))

avg_import_max_postprocessing= sum(import_max_postprocessing) / len(import_max_postprocessing)
BD_postprocessing=energy_demand_year/avg_import_max_postprocessing

charging_power = [P_ch[i].X for i in range(hours_in_year)]
discharging_power = [P_disch[i].X for i in range(hours_in_year)]
battery_energy = [E_b[i].X for i in range(hours_in_year + 1)]
selected_electricity_cost = [c_imp_low if high_usage.X <= 3500 else c_imp_high for _ in range(hours_in_year)]
#grid_usage_binary = [grid_usage[i].X for i in range(hours_in_year)]
#cumulative_hours_usage = np.cumsum(grid_usage_binary)
selected_battery_capacity = C_b.X  # Get the optimized battery capacity
cost_inv = C_b.X * UP_b / lifetime_b
cost_op_value = cost_op.getValue()
cost_peak_value = cost_peak.getValue()  

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
fig, axs = plt.subplots(2, 1, figsize=(12, 24))

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

plt.tight_layout()
plt.show()

# bar plots for cost

# Adjust the list of costs and labels
cost_categories = ['Investment Cost', 'Operational Cost', 'Peak Import Cost']
#cost_values = [float(cost_inv), float(cost_op), float(cost_peak.getValue())]  
cost_values = [cost_inv, cost_op_value, cost_peak_value]
 

# Bar Plot for Costs
plt.figure(figsize=(10, 6))
bar_positions = range(len(cost_categories))
plt.bar(bar_positions, cost_values, color=['skyblue', 'lightgreen', 'salmon'])
#plt.title('Investment vs Operational vs Peak Demand Costs')
plt.ylabel('Cost in CHF/y')
plt.xticks(bar_positions, cost_categories)

# Adding text labels for each bar with formatted cost values
for i, value in enumerate(cost_values):
    plt.text(i, value + 0.05 * max(cost_values), f'{value:.2f}', ha='center')

plt.tight_layout()
plt.show()

fig, axs = plt.subplots(3, 1, figsize=(12, 18))

# Plot 1 - Demand versus Time
axs[0].plot(energy_demand, label='Energy Demand (kWh)', color='blue')
axs[0].set_xlabel('Time (hours)')
axs[0].set_ylabel('Energy (kWh)')
axs[0].set_title('Energy Demand Over Time')
axs[0].legend()

# Plot 2 - Imported Power versus Time
axs[1].plot(import_power, label='Imported Power (kW)', color='green')
axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('Power (kW)')
axs[1].set_title('Imported Power Over Time')
axs[1].legend()
# Annotations for Plot 2
axs[1].annotate(f'Energy Demanded: {energy_demand_year:.2f} kWh/y', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'Imported Energy: {energy_imported_year:.2f} kWh/y', xy=(0.5, 0.8), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'Averaged Demanded Peak Power: {average_max_demand:.2f} kW/month', xy=(0.5, 0.7), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'Averaged Imported Peak Power - optimization: {avg_import_max_opt:.2f} kW/month', xy=(0.5, 0.6), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'Averaged Imported Peak Power - post processing: {avg_import_max_postprocessing:.2f} kW/month', xy=(0.5, 0.5), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'BD - opt: {BD_opt:.2f} hours', xy=(0.5, 0.4), xycoords='axes fraction', fontsize=13)
axs[1].annotate(f'BD - post prodcessing: {BD_postprocessing:.2f} hours', xy=(0.5, 0.3), xycoords='axes fraction', fontsize=13)

# Plot 3 - Energy Stored in the Battery versus Time
axs[2].plot(battery_energy[:-1], label='Battery Energy (kWh)', color='cyan')  # exclude last point for accurate indexing
axs[2].set_xlabel('Time (hours)')
axs[2].set_ylabel('Energy (kWh)')
axs[2].set_title('Battery Energy Over Time')
axs[2].legend()

plt.tight_layout()
plt.show()
