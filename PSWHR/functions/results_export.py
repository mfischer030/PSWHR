# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:03:05 2024

@author: fism
"""

# results_export.py

import os
import pandas as pd
from datetime import datetime
from gurobipy import GRB

def export_optimization_results(variable_names, results, results_directory):
    """
    Export optimization results to an Excel file.

    Parameters:
    - variable_names: List of variable names.
    - variable_values: Dictionary containing values of decision variables.
    - results_directory: Directory to save the Excel file.
    """
    # Extract optimization status
    optimization_status = results.get('status', None)

    # Check if the optimization was successful
    if optimization_status == GRB.OPTIMAL:
        # Create a DataFrame with the results
        results = pd.DataFrame({name: results[name] for name in variable_names})

        # Check if the results directory exists, create it if not
        if not os.path.exists(results_directory):
            os.makedirs(results_directory)

        # Get the current date and time
        current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Construct the Excel file name with date and hour
        excel_file_name = f'pv_sizingprob_results_{current_datetime}.xlsx'

        # Combine the directory and file name
        excel_file_path = os.path.join(results_directory, excel_file_name)

        # Save the DataFrame to an Excel file
        results.to_excel(excel_file_path, index=False)

        print(f"Results saved successfully at: {excel_file_path}")
    else:
        print("Export to excel failed")
