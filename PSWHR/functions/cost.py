# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:14:39 2024

@author: fism

"""
import pandas as pd
import numpy as np

def totalAnnualCost(system_sizes, energy_tariff,
                    UP, maintenance, life, 
                    P_imp, P_max_imp, P_exp, P_th_LT, P_th_HT,
                    cost_imp_el, cost_exp_el, cost_export_heatLT, cost_export_heatHT,
                    df_input, nHours, timeline_choice):
    
    # Capital recovery factor
    # CRF = (annual_interest_rate * (1 + annual_interest_rate)**project_lifetime) / ((1 + annual_interest_rate)**project_lifetime - 1) # not used

    # Installation and maintenance costs in [€/y]------------------------------
    cost_inst  = 0  # Initialize cost_inst outside the loop
    cost_maint = 0
    
    cost_inst  = sum((size * UP[component]) / life[component] for component, size in system_sizes.items())
    cost_maint = sum(size * UP[component] * maintenance[component] for component, size in system_sizes.items())
    
    # Operation costs in [€/y]-------------------------------------------------
    def electricity_prices(energy_tariff, P_imp, P_max_imp, P_exp, df_input, timeline_choice, vat_included=False):
        
        # Calculate the mean operating time (mittlere Benutzungsdauer)
        if timeline_choice == 'week':
            multiplier = 52                               # Assuming 52 weeks in a year
        elif timeline_choice == 'month':
            multiplier = 12                               # 12 months in a year
        else:
            multiplier = 1                                # Default to yearly calculation with no multiplication needed
        
        # # Calculate the mean of the monthly maximum values
        # mean_monthly_max = 0 
        # for month in P_max_imp: 
        #     mean_monthly_max += P_max_imp[month] / len(P_max_imp)
        
        # print(mean_monthly_max)
        # avr_operating_time = P_imp.sum() / mean_monthly_max
        
        avr_operating_time = 3600
        
        # print(f"The monthly maximum values is: {sum(P_imp)}")
        # print(f"The mean of the monthly maximum values is: {mean_monthly_max}")
        # print(f"The average operating time is: {avr_operating_time}")
        
        # Base fees in CHF for choosen timeline (week,month,year)
        base_fee_incl_vat = (42.16 + 616.17) / multiplier
        base_fee_excl_vat = (39.00 + 570.00) / multiplier

        # Tariff-specific prices in Rappen per kWh
        tariffs = {
            "Green": {"peak_excl_vat": 14.92, "peak_incl_vat": 16.13, "off_peak_excl_vat": 12.12, "off_peak_incl_vat": 13.10},
            "Blue":  {"peak_excl_vat": 12.42, "peak_incl_vat": 13.43, "off_peak_excl_vat": 9.62, "off_peak_incl_vat": 10.40},
            "Grey":  {"peak_excl_vat": 11.42, "peak_incl_vat": 12.35, "off_peak_excl_vat": 8.62, "off_peak_incl_vat": 9.32},
        }
        
        if avr_operating_time > 3500: 
            grid_usage_perkWh = {
                "Leistungstarif": {"excl_vat": 15.99, "incl_vat": 17.29},
                "Hochtarif":   {"excl_vat": 4.03, "incl_vat": 4.36},
                "Niedertarif": {"excl_vat": 2.01, "incl_vat": 2.17},
                }
        else:
            grid_usage_perkWh = { 
                "Leistungstarif": {"excl_vat": 7.55, "incl_vat": 8.16},
                "Hochtarif":   {"excl_vat": 7.47, "incl_vat": 8.08},
                "Niedertarif": {"excl_vat": 3.74, "incl_vat": 4.04},
                }
    
        # SystemServices_Swissgrid + Stromreserve + Gesetzliche_Foerderabgabe + Abgabe_an_die_Gemeinde
        additionnal_grid_fees = {"excl_vat": 0.75 + 1.20 + 2.30 + 1.50, 
                                "incl_vat" : 0.81 + 1.30 + 2.49 + 1.62
                                }

        if energy_tariff not in tariffs:
            raise ValueError("Invalid tariff. Choose Green, Blue, or Grey.")

        energy_tariff_choice = tariffs[energy_tariff]
        
        # Given prices for exported energy in Rp/kWh
        price_export_energy = [13.07, 7.73, 7.24, 8.66]
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
        price_export_energy = np.mean(price_export_energy)
        
        # Peak hours definition
        peak_time_start = 7
        peak_time_end = 21

        # Initializing costs
        cost_elec_imp = 0
        cost_elec_exp = 0 
        cost_elec     = base_fee_incl_vat if vat_included else base_fee_excl_vat
        
        for i, row in df_input.iterrows():
            hour = row['HR']
            if peak_time_start <= hour < peak_time_end:
                price_energy_per_kWh = energy_tariff_choice["peak_incl_vat" if vat_included else "peak_excl_vat"]
                price_grid_per_kWh = grid_usage_perkWh["Hochtarif"]["incl_vat" if vat_included else "excl_vat"]
                price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"] 
            else:
                price_energy_per_kWh = energy_tariff_choice["off_peak_incl_vat" if vat_included else "off_peak_excl_vat"]
                price_grid_per_kWh = grid_usage_perkWh["Niedertarif"]["incl_vat" if vat_included else "excl_vat"]
                price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"]
            
            cost_elec_imp += P_imp[i] / 1000 * price_energy_per_kWh / 100       # in € assuming 1CHF = 1€
            cost_elec_exp += P_exp[i] / 1000 * price_export_energy / 100        # in €
        
        # Calculating the grid fees associated with imports (Netzentgelte)
        sum_P_imp = 0
        sum_P_imp += sum(P_imp[t] for t in range(nHours))
        cost_grid_usage = sum_P_imp / 1000 * price_grid_per_kWh / 100 
        
        for month in P_max_imp:
            cost_grid_usage += P_max_imp[month] / 1000 * grid_usage_perkWh["Leistungstarif"]["incl_vat" if vat_included else "excl_vat"]
        
        # Total cost of electricity imports and exports + grid fees
        cost_elec += (cost_elec_imp - cost_elec_exp) + cost_grid_usage

        return cost_elec, cost_elec_imp, cost_elec_exp, cost_grid_usage
    
    cost_elec, cost_elec_imp, cost_elec_exp, cost_grid_usage = electricity_prices(energy_tariff, P_imp, P_max_imp, P_exp, df_input, timeline_choice, vat_included=False)
    
    cost_WHR = sum((P_th_LT[t]/1000 * cost_export_heatLT + P_th_HT[t]/1000 * cost_export_heatHT) for t in range(nHours))
    
    cost_op  = cost_elec - cost_WHR

    return cost_inst, cost_elec_imp, cost_elec_exp, cost_grid_usage, cost_elec, cost_op, cost_maint, cost_WHR

# Code from Roxanne------------------------------------------------------------

# cost without WHR
# cost_inst = (P_PV_peak*UP_PV/life_PV+P_e_nom*UP_e/life_e+C_b/3600*UP_b/life_b)/1000 # k€
# cost_op   = sum(P_imp.*cost_el)/1000 # k€

# costs with WHR
# cost_inst    = (P_PV_peak * UP_PV/life_PV + P_e_nom * UP_e/life_e + C_b/3600 * UP_b/life_b + S_HP * UP_HP/life_HP + 43872*(S_c/1000)^(0.5861)/life_c + UP_storage * mass_H2_day_obj/life_storage)/1000; # [kEUR]
# cost_op      = sum(P_imp.*cost_el)/1000 - sum(P_th_HT.*cost_export_heatHT)/1000 - sum(P_exp.*cost_export_el)/1000 # [kEUR]
# cost_startup = sum(startup*cost_startup_ely)/1000

# costs with more WHR and maintenance
# cost_inst    = (P_PV_peak * UP_PV * ann_PV + P_e_nom * UP_e * ann_e + C_b/3600 * UP_b * ann_b + S_HP * UP_HP * ann_HP + (A_HEX*UP_HEX + Fixed_HEX) * ann_HEX + 43872*(S_c)^(0.5861) * ann_c + UP_storage * mass_H2_day_obj * ann_storage + UP_disp)/1000; # [kEUR/y] add more components minutillo, dispenser and refrigeration, size known by size storage
# cost_inst    = (P_PV_peak * UP_PV * ann_PV + cost_e * ann_e + C_b/3600 * UP_b * ann_b + S_HP * UP_HP * ann_HP + (A_HEX*UP_HEX + Fixed_HEX) * ann_HEX + cost_c * ann_c + UP_storage * mass_H2_day_obj * ann_storage + UP_disp + P_refr * UP_refr * ann_refr )/1000; # [kEUR/y] add more components minutillo, dispenser and refrigeration
# cost_op      = sum([P_imp[i] * cost_el[i] for i in range(len(P_imp))]) / 1000 - sum([P_th_HT[i] * cost_export_heatHT[i] for i in range(len(P_th_HT))]) / 1000 - sum([P_th_LT[i] * cost_export_heatLT[i] for i in range(len(P_th_LT))]) / 1000 - sum([P_exp[i] * cost_export_el[i] for i in range(len(P_exp))]) / 1000                         # [kEUR/y]
# cost_maint   = (maint_PV * P_PV_peak * UP_PV + maint_e * cost_e + maint_b * C_b / 3600 * UP_b * ann_b + maint_HP * S_HP * UP_HP * ann_HP + maint_HEX * (A_HEX * UP_HEX + Fixed_HEX) * ann_HEX + maint_c * cost_c + maint_storage * UP_storage * mass_H2_day_obj * ann_storage + maint_disp * UP_disp + maint_refr * P_refr * UP_refr) / 1000  # [kEUR/y] mainly Poalo Gabrielli data
# cost_startup = sum([startup[i] * cost_startup_ely[i] for i in range(len(startup))]) / 1000
#------------------------------------------------------------------------------