# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:03:05 2024

@author: fism
"""

# results_export.py

from gurobipy import GRB
import os
import pandas as pd
from datetime import datetime
from openpyxl import Workbook

def export_optimization_results(variable_names, results, results_directory, all_costs, system_sizes):
    """
    Export optimization results along with all costs and system sizes to an Excel file.

    Parameters:
    - variable_names: List of variable names.
    - results: Dictionary containing the optimization status and values of decision variables.
    - results_directory: Directory to save the Excel file.
    - all_costs: Dictionary with all cost details.
    - system_sizes: Dictionary with system size details.
    """
    # Extract optimization status
    optimization_status = results.get('status', None)

    # Check if the optimization was successful
    if optimization_status == GRB.OPTIMAL:
        # Initialize DataFrame
        results_df = pd.DataFrame()

        # Flatten and convert the results dictionary to DataFrame
        for name in variable_names:
            if name in results:
                data = results[name]
                if isinstance(data, list):  # If the data is a list, expand it into rows
                    for i, value in enumerate(data):
                        results_df.at[i, name] = value
                else:
                    results_df[name] = [data]  # If not a list, just insert the data
        # Flatten and convert the results dictionary to DataFrame
        def flatten_dict(d, parent_key='', sep='_'):
            items = []
            for k, v in d.items():
                new_key = f"{parent_key}{sep}{k}" if parent_key else k
                if isinstance(v, dict):
                    items.extend(flatten_dict(v, new_key, sep=sep).items())
                else:
                    items.append((new_key, v))
            return dict(items)

        # Assume the flatten_dict function has been defined earlier
        all_costs_df = pd.DataFrame([flatten_dict(all_costs)])
        system_sizes_df = pd.DataFrame([flatten_dict(system_sizes)])

        # Check if the results directory exists, create it if not
        if not os.path.exists(results_directory):
            os.makedirs(results_directory)

        # Get the current date and time
        current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Construct the Excel file name with date and time
        excel_file_name = f'optimization_results_{current_datetime}.xlsx'
        excel_file_path = os.path.join(results_directory, excel_file_name)

        # Save all DataFrames to the same Excel file with different sheets
        with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
            results_df.to_excel(writer, sheet_name='Optimization Results', index=False)
            all_costs_df.to_excel(writer, sheet_name='Cost', index=False)
            system_sizes_df.to_excel(writer, sheet_name='Sizes', index=False)

        print(f"Results saved successfully at: {excel_file_path}")
    else:
        print("Export to excel failed: Optimization was not optimal or status was missing.")