# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:13:11 2024

@author: fism

PWA approximation of PEM Electrolyser and Fuel Cell efficiency
"""

import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
# PV
#------------------------------------------------------------------------------
"""Sources: De Soto W et al. (2006); Sun V et al. (2020); Dubey et al. (2013)"""

def pv_efficiency(irradiance, T_amb):
    eta_PV_ref  = 0.15  # Consider changing back to 0.21 
    T_PV_ref    = 25    # Reference temperature of the pv cell in [°C]
    beta_PV_ref = 0.004 # in [K]
    gamma_PV    = 0.12
    NOCT        = 45
    
    # Calculate T_cell based on ambient Temperature T_amb (also found in 2015_Guinot et al.)
    T_cell = T_amb + (NOCT - 20) * irradiance / 800

    # Initialize eff_cell as a zeros array of the same length as irradiance
    eta_cell = np.zeros(len(irradiance))

    # Loop through each irradiance value
    for i in range(len(irradiance)): 
        if irradiance[i] == 0:
            eta_cell[i] = 0
        else:
            eta_cell[i] = eta_PV_ref * (1 - beta_PV_ref * (T_cell[i] - T_PV_ref) + gamma_PV * np.log10(irradiance[i]))
    
    return eta_cell

#------------------------------------------------------------------------------
# ELECTROLYSER
#------------------------------------------------------------------------------
"""
Using Paolo Gabrielli's script, following is found to describe the
efficiency of the electrolyser, using the HHV of hydrogen (iso LHV in
previous script)
for electrolyser of size 100kW, normalized
"""

def pwa_eta_ELY(N_bp):
    if N_bp == 2:
        x_bp_val = np.array([0.069264723, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.584970932])
        mm_elec  = 0.544695364
        qq_elec  = -0.009873526  

    elif N_bp == 3:
        x_bp_val = np.array([0.069264723, 0.56492096, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.312682424, 0.584970932])
        mm_elec  = np.array([0.574647821, 0.51653222])
        qq_elec  = np.array([-0.011948175, 0.020882546]) 

    elif N_bp == 4:
        x_bp_val = np.array([0.069264723, 0.378789088, 0.735922253, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.210531814, 0.403454818,	0.584970932])
        mm_elec = np.array([0.59018671,	0.540199071, 0.509667852])
        qq_elec = np.array([-0.013024472, 0.005910301,	0.028378904]) 
    
    elif N_bp == 8:
        x_bp_val = np.array([0.069264723, 0.200398341, 0.338414769, 0.481335201,
                             0.628463327, 0.77943386, 0.934020342, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.10825737, 0.187821926, 0.267288394,
                             0.34672396, 0.426145867, 0.505560544, 0.584970932])
        mm_elec = np.array([0.613135855, 0.57648612, 0.556018952, 0.539907411,
                            0.526075556, 0.513723294, 0.502445365])
        qq_elec = np.array([-0.014614038, -0.007269492, -0.0003431, 0.007411952,
                            0.016104766, 0.025732536, 0.036266352])
    
    # Coefficient from Roxanne for S_ELY=400kW
    
    # if N_bp == 2:
    #     x_bp_val = np.array([77.6555, 428.5104])
    #     y_bp_val = np.array([56.5484, 276.4963])
    #     mm_elec = 0.6269
    #     qq_elec = 7.8669 / 400  # Normalizing the intersection with y-axis, see Paolo Gabrielli

    # elif N_bp == 3:
    #     x_bp_val = np.array([77.6555, 250.0070, 428.5104])
    #     y_bp_val = np.array([56.5484, 169.2480, 276.4963])
    #     mm_elec = np.array([0.6539, 0.6008])
    #     qq_elec = np.array([5.7700, 19.0390]) / 400

    # elif N_bp == 4:
    #     x_bp_val = np.array([77.6555, 190.6689, 302.2049, 428.5104])
    #     y_bp_val = np.array([56.5484, 131.7027, 201.4249, 276.4963])
    #     mm_elec = np.array([0.6650, 0.6251, 0.5944])
    #     qq_elec = np.array([4.9073, 12.5138, 21.8054]) / 400
       
        # poly_eff_380kW = np.array([-0.0001722, 0.7099, 2.0])    # For S_ELY 380 kW
        # poly_eff_510kW = np.array([-0.0001305, 0.7099, 2.95])     # For S_ELY 510 kW
        # poly_eff_187kW = np.array([-0.0003559, 0.7099, 1.083])  # For S_ELY 187 kW
    
    # elif N_bp == 5:
    #     x_bp_val = np.array([77.6555, 157.5608, 250.0070, 337.6598, 428.5104])
    #     y_bp_val = np.array([56.5484, 110.2433, 169.2480, 222.8748, 276.4963])
    #     mm_elec = np.array([0.6720, 0.6383, 0.6118, 0.5902])
    #     qq_elec = np.array([4.3654, 9.6786, 16.2915, 23.5824]) / 400
    
    # elif N_bp == 7:
    #     x_bp_val = np.array([77.6555, 133.1374, 190.6689, 250.0070, 302.2049, 364.5867, 428.5104])
    #     mm_elec = np.array([0.6776, 0.6528, 0.6327, 0.6164, 0.6017, 0.5872])
    #     qq_elec = np.array([3.9268, 7.2290, 11.0598, 15.1335, 19.5806, 24.8839]) / 400
    
    # elif N_bp == 10:
    #     x_bp_val = np.array([77.6555, 117.0607, 149.3797, 190.6689, 224.3651, 267.2727, 302.2049, 346.6037, 382.6954, 428.5104])
    #     y_bp_val = np.array([56.5484, 83.4089, 104.8776, 131.7027, 153.1581, 179.9740, 201.4249, 228.2371, 249.6860, 276.4963])
    #     mm_elec = np.array([0.6816, 0.6643, 0.6497, 0.6367, 0.6250, 0.6141, 0.6039, 0.5943, 0.5852])
    #     qq_elec = np.array([3.6146, 5.6488, 7.8273, 10.2979, 12.9370, 15.8493, 18.9253, 22.2545, 25.7377]) / 400
    
    # elif N_bp == 12:
    #     x_bp_val = np.array([77.6555, 109.0870, 133.1374, 165.7809, 199.0390, 232.8778, 258.6230, 293.4224, 328.7478, 364.5867, 391.7966, 428.5104])
    #     mm_elec = np.array([0.6838, 0.6696, 0.6575, 0.6452, 0.6340, 0.6250, 0.6164, 0.6072, 0.5985, 0.5912, 0.5842])
    #     qq_elec = np.array([3.4513, 4.9922, 6.6029, 8.6465, 10.8708, 12.9838, 15.1892, 17.8923, 20.7603, 23.4182, 26.1607]) / 400

    return x_bp_val, y_bp_val, mm_elec, qq_elec

#------------------------------------------------------------------------------
# FUEL CELL
#------------------------------------------------------------------------------

param_struct = {
    'P_max': 80,
    'P_min': 30,
    'fc_i_breakpoint': 123.054000168826,       # mA
    'fc_p_breakpoint': 47.9673651830867,       # kW
    'fc_slope_high_current': 3.31471893412921,
    'fc_slope_low_current': 2.56536917754678,
    'fc_i_to_h2': 3.52555108268819,    # vdot calculation explained in the comment
    'fc_dc_to_ac': 0.956183048244058   # relation between dc input and ac power delivered
}

# Parameters provided
P_bp     = 47.97 * 1000         # Breakpoint power in Watts (converted from kW)
i_bp     = 123.05               # Breakpoint current in A
P_min_fc = 30 * 1000            # Minimum power in Watts (converted from kW)
P_max_fc = 80 * 1000            # Maximum power in Watts (converted from kW)
s1       = 2.56                 # Slope below breakpoint
s2       = 3.31                 # Slope above breakpoint
c        = 0.21                 # Volume flow conversion factor
d        = 0.96                 # Efficiency coefficient


# fc_i_to_h2 calculation
# vdot = #cells * M / (zF * rho) * 3600s/h * 1000 Nl/Nm3 * 1 / (60min/h) * i
#      = #cells * 0.007 Nl/(A*min) * i 

 # 0.007 vs 3.525?? -> 3.525 / 0.007 = 503.6502

# param_struct['fc_i_to_o2'] = 1.85306635062262

#------------------------------------------------------------------------------
# COMPRESSOR
#------------------------------------------------------------------------------
"""
Sources: thesis (Josien de Koning and) Timo Laaksonlaita
"""

def compressor_curve():
    mm_c = np.array([4.8756, 6.5348, 2.5470])
    med_c = 4.078  # Midway H2 flow compressor [kg/h]
    return mm_c, med_c

