To get the code running please adjust the following links in the code once you downloaded the pvsizingprob folder:
1. System path: Please adjust the following code in PV_SIZINGPROB (Line 31):
     - here should be the path to the directory where the functions import_data, plotting_module and results_export should be!

"""
path to the module 'import_data'

sys.path.append(r'C:\Users\...\02_modeling_and_optimization\functions') 
"""

2. FUNCTION > IMPORT_DATA.PY
    - Line 118: update link to: Inputs_BoundaryLoads_week.xlsx
    - Line 119: update link to: demand_data_weekly.xlsx
    - Line 128: update link to: Rubigen_2019-2022_irradiance_Hourly.xlsx
    - Line 129: update link to: demand_data_yearly.xlsx
3. PV_SIZINGPROB.PY
    - Update path where you want to save the results > I already created a directory named 'results' for this purpose in the pvsizingprob folder.
      
    """
    Define the path to the results directory
   
    results_directory = r'C:\Users\...\pvsizingprob\results'
    """
