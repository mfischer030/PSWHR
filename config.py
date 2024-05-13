# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 14:44:23 2024

Configuration script - create a new user by following the given structure. 
The following 4 Links need to be addressed in the paths_configuration function:

input_path  = path to irradiance data (Rubigen_2019-2022_irradiance_Hourly.xlsx)
demand_path = path to demand profiles (200923_Lastang 2015-2020_Rubigen.xlsx)
function_path = path to the folder with all the fucntions (...\PSWHR\functions)
export_path = path to the directory where to export results 

@author: fism
"""

def paths_configuration(user):
    if user == 'maxime':
        input_path  = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\PSWHR\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx'
        demand_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\PSWHR\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\PSWHR\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\PSWHR\PSWHR\functions'
        export_path = r'C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\results'
    
    elif user == 'maxime_EMPA_WS':
        input_path  = r'C:\Users\fism\Desktop\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx'
        demand_path = r'C:\Users\fism\Desktop\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path = r'C:\Users\fism\Desktop\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\fism\Desktop\PSWHR\functions'
        export_path = r'C:\Users\fism\Desktop\PSWHR\results'
    
    elif user == 'christian':
        input_path  = r'C:\Users\peter_c\Desktop\test\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
        demand_path = r'C:\Users\peter_c\Desktop\test\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path   = r'C:\Users\peter_c\Desktop\test\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\peter_c\Desktop\test\PSWHR\functions'
        export_path = r'C:\Users\peter_c\Desktop\test\results'
    
    # elif user == 'gabriele':
        input_path  = r'C:\Users\...\PSWHR\input_data\Rubigen_2019-2022_irradiance_Hourly.xlsx',
        demand_path = r'C:\Users\...\PSWHR\input_data\200923_Lastang 2015-2020_Rubigen.xlsx'
        heat_path   = r'C:\Users\...\PSWHR\input_data\Heat_Inputs_year_2022.xlsx'
        function_path = r'C:\Users\...\PSWHR\functions'
        export_path = r'C:\Users\...\results'
    
    
    return input_path, demand_path, heat_path, function_path, export_path
