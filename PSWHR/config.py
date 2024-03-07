# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 14:44:23 2024

Configuration script

input_path  = path to irradiance data (Rubigen_2019-2022_irradiance_Hourly.xlsx)
demand_path = path to demand profiles (200923_Lastang 2015-2020_Rubigen.xlsx)
function_path = path to the folder with all the fucntions
export_path = path to the directory where to export results 

@author: fism
"""

def paths_configuration(user):
    if user == 'maxime':
        input_path  = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\Rubigen_2019-2022_irradiance_Hourly.xlsx'
        demand_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\inputData\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\functions'
        export_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\results'
    
    elif user == 'christian':
        input_path  = r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
        demand_path = r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path   = r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\peter_c\Desktop\test\PSWHR\PSWHR\functions'
        export_path = r'C:\Users\peter_c\Desktop\test\results'
    
    # elif user == 'gabriele':
        input_path  = r'C:\Users\...\PSWHR\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
        demand_path = r'C:\Users\...\PSWHR\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path   = r'C:\Users\...\PSWHR\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\...\PSWHR\PSWHR\functions'
        export_path = r'C:\Users\...\results'
    
    
    return input_path, demand_path, heat_path, function_path, export_path
