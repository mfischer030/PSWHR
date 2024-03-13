# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 12:15:40 2023

@author: fism
"""
# Function calculates the annual electricity cost for industrial 
# customers with a yearly consumption over 50,000 kWh, considering different 
# tariff options (Green, Blue, Grey) and distinguishing between peak 
# and off-peak hours. 
# The user can specify whether the price includes Value Added Tax (VAT) or not.
# URL:https://www.bkw.ch/de/energie/gesetzliche-publikationen/tarife-tarifanpassungen
#------------------------------------------------------------------------------
import pandas as pd
import numpy as np

input_path  = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\Rubigen_2019-2022_irradiance_Hourly.xlsx'
demand_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\200923_Lastang 2015-2020_Rubigen.xlsx'

results = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\results\pv_sizingprob_results_2024-03-04_16-14-05.xlsx'   # year
#results = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\results\pv_sizingprob_results_2024-02-14_18-03-17.xlsx'  # month


df_input = pd.read_excel(input_path, sheet_name="2019")              # year
#df_input = pd.read_excel(input_path, sheet_name="2019_month_win")   # month
df_results = pd.read_excel(results)

p_imp = df_results['P_imp'].values
p_exp = df_results['P_exp'].values

def electricity_prices(energy_tariff, p_imp, p_exp, df_input, timeline_choice, vat_included=True):
    
    # Create DataFrame for p_imp with the MO list indicating the month
    df_imp = pd.DataFrame({'Month': df_input['MO'], 'p_imp': p_imp})

    # Group by month and find the maximum for each month for p_imp and p_exp
    monthly_max_imp = df_imp.groupby('Month')['p_imp'].max()
    
    # Calculate the mean of the monthly maximum values
    mean_monthly_max = monthly_max_imp.mean() 
    
    # Calculate the mean operating time (mittlere Benutzungsdauer)
    if timeline_choice == 'week':
        multiplier = 52                               # Assuming 52 weeks in a year
    elif timeline_choice == 'month':
        multiplier = 12                               # 12 months in a year
    else:
        multiplier = 1                                # Default to yearly calculation with no multiplication needed
    
    avr_operating_time = (sum(p_imp) / mean_monthly_max) * multiplier
    
    print(f"The monthly maximum values is: {sum(p_imp)}")
    print(f"The mean of the monthly maximum values is: {mean_monthly_max}")
    print(f"The average operating time is: {avr_operating_time}")
    
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
            #"Leistungstarif": {"excl_vat": 15.99, "incl_vat": 17.29},
            "Hochtarif":   {"excl_vat": 4.03, "incl_vat": 4.36},
            "Niedertarif": {"excl_vat": 2.01, "incl_vat": 2.17},
            }
    else:
        grid_usage_perkWh = { 
            #"Leistungstarif": {"excl_vat": 7.55, "incl_vat": 8.16},
            "Hochtarif":   {"excl_vat": 7.47, "incl_vat": 8.08},
            "Niedertarif": {"excl_vat": 3.74, "incl_vat": 4.04},
            }
    
    # additionnal_grid_fees = {
    #     "SystemServices_Swissgrid" : {"excl_MWSt": 0.75, "incl_vat": 0.81},
    #     "Stromreserve"             : {"excl_MWSt": 1.20, "incl_vat": 1.30},
    #     "Gesetzliche_Foerderabgabe": {"excl_MWSt": 2.30, "incl_vat": 2.49},
    #     "Abgabe_an_die_Gemeinde"   : {"excl_MWSt": 1.50, "incl_vat": 1.62},
    #     }
    additionnal_grid_fees = {
        "excl_vat": 0.75 + 1.20 + 2.30 + 1.50, 
        "incl_vat" : 0.81 + 1.30 + 2.49 + 1.62
        }
    
    # if grid_tariff == "Hochtarif":
    #     price_grid_per_kWh = grid_usage_perkWh["Hochtarif"]["incl_vat" if vat_included else "excl_vat"]
    #     price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"] 
        
    # elif grid_tariff == "Niedertarif":
    #     price_grid_per_kWh = grid_usage_perkWh["Niedertarif"]["incl_vat" if vat_included else "excl_vat"]
    #     price_grid_per_kWh += additionnal_grid_fees["incl_vat" if vat_included else "excl_vat"]
        

    if energy_tariff not in tariffs:
        raise ValueError("Invalid tariff. Choose Green, Blue, or Grey.")

    energy_tariff_choice = tariffs[energy_tariff]
    
    # Given prices for exported energy
    price_exp_energy = [13.07, 7.73, 7.24, 8.66]

    # Calculate the mean of these prices using numpy
    price_exp_energy = np.mean(price_exp_energy)
    
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
        
        cost_elec_imp += p_imp[i] / 1000 * price_energy_per_kWh / 100           # in € assuming 1CHF = 1€
        cost_elec_exp += p_exp[i] / 1000 * price_exp_energy / 100  # in €
    
    cost_grid_usage = (sum(p_imp) / 1000) * price_grid_per_kWh / 100
    
    cost_elec += (cost_elec_imp - cost_elec_exp) + cost_grid_usage

    return cost_elec, cost_elec_imp, cost_elec_exp, cost_grid_usage


energy_tariff = "Blue"
timeline_choice = 'year'


cost_elec = electricity_prices(energy_tariff, p_imp, p_exp, df_input, timeline_choice, vat_included=True)
print(cost_elec)
