# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:14:39 2024

@author: fism

"""
import pandas as pd
import numpy as np

def totalAnnualCost(system_sizes, energy_tariff, discountRate,
                    UP, maintenance, life, 
                    P_imp, P_max_imp, P_exp, P_th_LT, P_th_HT,
                    cost_imp_el, cost_exp_el, cost_export_heatLT, cost_export_heatHT,
                    m, high_usage,
                    df_input, nHours, timeline_choice,
                    cost_ELY, #cost_C
                    ):
    

    # Installation and maintenance costs in [€/y]------------------------------
    cost_inst  = 0  # Initialize cost_inst outside the loop
    cost_maint = 0
    
    r = discountRate
    cost_inst_PV  = system_sizes['PV'] * UP['PV'] *  r / (1 - (1 / ((1 + r) ** life['PV'])))
    cost_maint_PV = system_sizes['PV'] * UP['PV'] *  r / (1 - (1 / ((1 + r) ** life['PV']))) * maintenance['PV']
    
    cost_inst_BAT  = system_sizes['BAT'] * UP['BAT'] *  r / (1 - (1 / ((1 + r) ** life['BAT'])))
    cost_maint_BAT = system_sizes['BAT'] * UP['BAT'] *  r / (1 - (1 / ((1 + r) ** life['BAT']))) * maintenance['BAT']
    
    cost_inst_ELY  = cost_ELY *  r / (1 - (1 / ((1 + r) ** life['ELY'])))
    cost_maint_ELY = cost_ELY *  r / (1 - (1 / ((1 + r) ** life['ELY']))) * maintenance['ELY']
    
    cost_inst_C  = system_sizes['HESS']['C'] * UP['C'] *  r / (1 - (1 / ((1 + r) ** life['C'])))
    cost_maint_C = system_sizes['HESS']['C'] * UP['C'] *  r / (1 - (1 / ((1 + r) ** life['C']))) * maintenance['C']
    
    # cost_inst_C  = cost_C *  r / (1 - (1 / ((1 + r) ** life['C'])))
    # cost_maint_C = cost_C *  r / (1 - (1 / ((1 + r) ** life['C']))) * maintenance['C']
    
    cost_inst_TANK  = system_sizes['HESS']['TANK'] * UP['TANK'] *  r / (1 - (1 / ((1 + r) ** life['TANK'])))
    cost_maint_TANK = system_sizes['HESS']['TANK'] * UP['TANK'] *  r / (1 - (1 / ((1 + r) ** life['TANK']))) * maintenance['TANK']
    
    cost_inst_FC  = system_sizes['HESS']['FC'] * UP['FC'] *  r / (1 - (1 / ((1 + r) ** life['FC'])))
    cost_maint_FC = system_sizes['HESS']['FC'] * UP['FC'] *  r / (1 - (1 / ((1 + r) ** life['FC']))) * maintenance['FC']
    
    cost_inst_HEX  = system_sizes['WHR']['HEX'] * UP['HEX'] *  r / (1 - (1 / ((1 + r) ** life['HEX'])))
    cost_maint_HEX = system_sizes['WHR']['HEX'] * UP['HEX'] *  r / (1 - (1 / ((1 + r) ** life['HEX']))) * maintenance['HEX']
    
    cost_inst_HP  = system_sizes['WHR']['HP'] * UP['HP'] *  r / (1 - (1 / ((1 + r) ** life['HP'])))
    cost_maint_HP = system_sizes['WHR']['HP'] * UP['HP'] *  r / (1 - (1 / ((1 + r) ** life['HP']))) * maintenance['HP']
    
    cost_inst  = cost_inst_BAT + cost_inst_C + cost_inst_ELY + cost_inst_FC + cost_inst_HEX + cost_inst_HP + cost_inst_TANK + cost_inst_PV
    cost_maint = cost_maint_BAT + cost_maint_C + cost_maint_ELY + cost_maint_FC + cost_maint_HEX + cost_maint_HP + cost_maint_TANK + cost_maint_PV
    
    # cost_inst  = sum((size * UP[component]) / life[component] for component, size in system_sizes.items())
    # cost_maint = sum(size * UP[component] * maintenance[component] for component, size in system_sizes.items())
    
   #  def calculate_installation_cost_with_annuity(group, data, UP, life, annuityFactor):
   #      """Calculate installation costs using the annuity factor for amortization over component life."""
   #      if isinstance(data[group], dict):  # Handling sub-dictionaries like HESS and WHR
   #          return sum((size * UP[component] * annuityFactor[component]) for component, size in data[group].items())
   #      else:  # Handling single components like PV and BAT
   #          return (data[group] * UP[group] * annuityFactor[group])

   #  def calculate_maintenance_cost(group, data, UP, maintenance):
   #      """Calculate maintenance costs for a group of components or a single component."""
   #      if isinstance(data[group], dict):  # Handling sub-dictionaries like HESS and WHR
   #          return sum((size * UP[component] * maintenance[component]) for component, size in data[group].items())
   #      else:  # Handling single components like PV and BAT
   #          return (data[group] * UP[group] * maintenance[group])
    
   #  def detailed_costs_for_group_with_annuity(group_data, UP, life, maintenance, annuityFactor):
   #      """Calculate detailed installation and maintenance costs for components within a group using annuity factors."""
   #      return {
   #          component: {
   #              'installation': (size * UP[component] * annuityFactor[component]),
   #              'maintenance': (size * UP[component] * maintenance[component])
   #          }
   #          for component, size in group_data.items()
   #      }
    
   #  r = discountRate  # Annual Discount Rate
   #  annuityFactor = {
   #      component: r / (1 - (1 / ((1 + r) ** life[component])))
   #      for component in life
   # }
    
   #  # Calculate costs for each group, including detailed costs for HESS and WHR
   #  costs_with_annuity = {
   #      'PV': {
   #          'installation': calculate_installation_cost_with_annuity('PV', system_sizes, UP, life, annuityFactor),
   #          'maintenance': calculate_maintenance_cost('PV', system_sizes, UP, maintenance)
   #      },
   #      'BAT': {
   #          'installation': calculate_installation_cost_with_annuity('BAT', system_sizes, UP, life, annuityFactor),
   #          'maintenance': calculate_maintenance_cost('BAT', system_sizes, UP, maintenance)
   #      },
   #      'HESS': detailed_costs_for_group_with_annuity(system_sizes['HESS'], UP, life, maintenance, annuityFactor),
   #      'WHR': detailed_costs_for_group_with_annuity(system_sizes['WHR'], UP, life, maintenance, annuityFactor)
   #  }
    
   #  # Print the detailed costs dictionary
   #  print(costs_with_annuity)
    
   #  def calculate_total_costs(costs):
   #      """Calculate total installation and maintenance costs from a structured costs dictionary."""
   #      total_installation_cost = 0
   #      total_maintenance_cost = 0
   #      for group, details in costs.items():
   #          if isinstance(details, dict):
   #              # Check if direct group costs exist like PV and BAT
   #              if 'installation' in details and 'maintenance' in details:
   #                  total_installation_cost += details['installation']
   #                  total_maintenance_cost += details['maintenance']
   #              else:
   #                  # Sub-groups like HESS and WHR
   #                  for component, comp_details in details.items():
   #                      total_installation_cost += comp_details['installation']
   #                      total_maintenance_cost += comp_details['maintenance']
   #      return total_installation_cost, total_maintenance_cost

   #  # Calculate both total installation and maintenance costs
   #  cost_inst, cost_maint = calculate_total_costs(costs_with_annuity)
    
    print(f"Total Installation Cost: {cost_inst}")
    print(f"Total Maintenance Cost: {cost_maint}")
    
    # Operation costs in [€/y]-------------------------------------------------
    def electricity_cost(energy_tariff, P_imp, P_max_imp, P_exp, m, high_usage, df_input, timeline_choice, vat_included=True):
        
        # Base fees in CHF for choosen timeline (week,month,year)
        base_fee_incl_vat = 42.16 + 616.17 
        base_fee_excl_vat = 39.00 + 570.00 

        # Tariff-specific prices in Rappen per kWh
        tariffs = {
            "Green": {"peak_excl_vat": 14.92, "peak_incl_vat": 16.13, "off_peak_excl_vat": 12.12, "off_peak_incl_vat": 13.10},
            "Blue":  {"peak_excl_vat": 12.42, "peak_incl_vat": 13.43, "off_peak_excl_vat": 9.62, "off_peak_incl_vat": 10.40},
            "Grey":  {"peak_excl_vat": 11.42, "peak_incl_vat": 12.35, "off_peak_excl_vat": 8.62, "off_peak_incl_vat": 9.32},
        }
        
        if m.addConstr(high_usage >= 1, "BD>3500"): 
            grid_usage_perkWh = {
                "Leistungstarif": {"excl_vat": 15.99, "incl_vat": 17.29},
                "Hochtarif":   {"excl_vat": 4.03, "incl_vat": 4.36},
                "Niedertarif": {"excl_vat": 2.01, "incl_vat": 2.17},
                }
        elif m.addConstr(high_usage < 1, "BD<3500"):
            grid_usage_perkWh = { 
                "Leistungstarif": {"excl_vat": 7.55, "incl_vat": 8.16},
                "Hochtarif":   {"excl_vat": 7.47, "incl_vat": 8.08},
                "Niedertarif": {"excl_vat": 3.74, "incl_vat": 4.04},
                }
        
        # Use these to check BIG-M constraint: high values for BD>3500 should
        # result in the optimizer choosing prices for BD<3500
        
        # if m.addConstr(high_usage >= 1, "BD>3500"): 
        #     grid_usage_perkWh = {
        #         "Leistungstarif": {"excl_vat": 100, "incl_vat": 100},
        #         "Hochtarif":   {"excl_vat": 100, "incl_vat": 100},
        #         "Niedertarif": {"excl_vat": 100, "incl_vat": 100},
        #         }
        # elif m.addConstr(high_usage < 1, "BD<3500"):
        #     grid_usage_perkWh = { 
        #         "Leistungstarif": {"excl_vat": 7.55, "incl_vat": 8.16},
        #         "Hochtarif":   {"excl_vat": 7.47, "incl_vat": 8.08},
        #         "Niedertarif": {"excl_vat": 3.74, "incl_vat": 4.04},
        #         }
    
        # SystemServices_Swissgrid + Stromreserve + Gesetzliche_Foerderabgabe + Abgabe_an_die_Gemeinde
        additionnal_grid_fees = {"excl_vat": 0.75 + 1.20 + 2.30 + 1.50, 
                                "incl_vat" : 0.81 + 1.30 + 2.49 + 1.62
                                }

        if energy_tariff not in tariffs:
            raise ValueError("Invalid tariff. Choose Green, Blue, or Grey.")

        energy_tariff_choice = tariffs[energy_tariff]
        
        # Given prices for exported energy in Rp/kWh
        price_export_energy = [13.07, 7.73, 7.24, 8.66]
        # price_export_energy = 0
        
        # Function to determine quartal based on month
        def determine_quartal(month):
            if 1 <= month <= 3: 
                return 0  # First quartal
            elif 4 <= month <= 6:
                return 1  # Second quartal
            elif 7 <= month <= 9:
                return 2  # Third quartal
            else:
                return 3  # Fourth quartal

        # Map each month to its corresponding quartal
        df_input['Quartal'] = df_input['MO'].apply(determine_quartal)

        # Assign the price of exported electricity based on the quartal
        df_input['Price_Export_Energy'] = df_input['Quartal'].apply(lambda x: price_export_energy[x])

        # Calculate the mean of these prices using numpy
        price_export_energy = df_input['Price_Export_Energy']
        
        # Peak hours definition
        peak_time_start = 7
        peak_time_end = 21

        # Initializing costs
        cost_elec_imp = 0
        cost_elec_exp = 0 
        cost_grid_usage = 0
        cost_elec     = base_fee_incl_vat if vat_included else base_fee_excl_vat
        
        price_energy_per_kWh = []
        price_grid_per_kWh   = []
        
        for i, row in df_input.iterrows():
            hour = row['HR']
            if peak_time_start <= hour < peak_time_end:
                hourly_price_energy_per_kWh = energy_tariff_choice["peak_incl_vat" if vat_included else "peak_excl_vat"]
                hourly_price_grid_per_kWh = grid_usage_perkWh["Hochtarif"]["incl_vat" if vat_included else "excl_vat"]
                hourly_price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"] 
            else:
                hourly_price_energy_per_kWh = energy_tariff_choice["off_peak_incl_vat" if vat_included else "off_peak_excl_vat"]
                hourly_price_grid_per_kWh = grid_usage_perkWh["Niedertarif"]["incl_vat" if vat_included else "excl_vat"]
                hourly_price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"]
            
            # Create a list with prices for plotting
            price_energy_per_kWh.append(hourly_price_energy_per_kWh)
            price_grid_per_kWh.append(hourly_price_grid_per_kWh)
            
            # Determine the cost / revenues for importing and exporting electricity
            cost_elec_imp += P_imp[i] / 1000 * hourly_price_energy_per_kWh / 100       # in CHF
            cost_elec_exp += P_exp[i] / 1000 * price_export_energy[i] / 100            # in CHF
        
            # Calculating the grid fees associated with imports (Netzentgelte)
            # sum_P_imp = 0
            # sum_P_imp += sum(P_imp[t] for t in range(nHours))
            cost_grid_usage += P_imp[i] / 1000 * hourly_price_grid_per_kWh / 100 
            
        
        for month in P_max_imp:
            cost_grid_usage += P_max_imp[month] / 1000 * grid_usage_perkWh["Leistungstarif"]["incl_vat" if vat_included else "excl_vat"]
        
        # Total cost of electricity imports and exports + grid fees
        cost_elec += (cost_elec_imp - cost_elec_exp) + cost_grid_usage
        
        electricity_prices = {
            'Electricity prices [Rp./kWh]': price_energy_per_kWh,
            'Grid use cost [Rp./kWh]': price_grid_per_kWh,
            'Fees for monthly Import peaks [CHF/kW/Monat]':grid_usage_perkWh["Leistungstarif"],
            'Additional Fees [Rp/kWh] (already accounted in Grid use cost)': additionnal_grid_fees,
            'Elecricity export price [Rp.kWh]': price_export_energy,
            'Selected tariffs for energy prices': tariffs,
            'Selected tariffs for grid usage prices': grid_usage_perkWh
            }

        return cost_elec, cost_elec_imp, cost_elec_exp, cost_grid_usage, electricity_prices, 
    
    
    cost_elec, cost_elec_imp, cost_elec_exp, cost_grid_usage, electricity_prices = electricity_cost(energy_tariff, P_imp, P_max_imp, P_exp, m, high_usage, df_input, timeline_choice, vat_included=True)
    
    # Calculate the revenues of exported Heat, prices for heat are in [€/kWh]
    cost_WHR = sum((P_th_LT[t]/1000 * cost_export_heatLT + P_th_HT[t]/1000 * cost_export_heatHT) for t in range(nHours))
    
    # Calculate the total operating cost
    cost_op  = cost_elec - cost_WHR

    return cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, cost_maint, cost_WHR, electricity_prices, #costs_with_annuity

