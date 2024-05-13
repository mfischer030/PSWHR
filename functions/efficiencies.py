# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:13:11 2024

@author: fism

PWA approximation of PEM Electrolyser and Fuel Cell efficiency
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# PV
# ------------------------------------------------------------------------------
"""Sources: De Soto W et al. (2006); Sun V et al. (2020); Dubey et al. (2013)"""

def calculate_pv_efficiency(irradiance, T_amb):
    eta_PV_ref = 0.15  # Consider changing back to 0.21
    T_PV_ref = 25    # Reference temperature of the pv cell in [°C]
    beta_PV_ref = 0.004  # in [K]
    gamma_PV = 0.12
    NOCT = 45

    # Calculate T_cell based on ambient Temperature T_amb (also found in 2015_Guinot et al.)
    T_cell = T_amb + (NOCT - 20) * irradiance / 800

    # Initialize eff_cell as a zeros array of the same length as irradiance
    eta_cell = np.zeros(len(irradiance))

    # Loop through each irradiance value
    for i in range(len(irradiance)):
        if irradiance[i] == 0:
            eta_cell[i] = 0
        else:
            eta_cell[i] = eta_PV_ref * (1 - beta_PV_ref * (
                T_cell[i] - T_PV_ref) + gamma_PV * np.log10(irradiance[i]))

    return eta_cell

# ------------------------------------------------------------------------------
# ELECTROLYSER
# ------------------------------------------------------------------------------
"""
Using Paolo Gabrielli's script, following is found to describe the
efficiency of the electrolyser, using the HHV of hydrogen (iso LHV in
previous script)
for electrolyser of size 100kW, normalized
"""

def pwa_ELY(N_bp):
    if N_bp == 2:
        x_bp_val = np.array([0.069264723, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.584970932])
        mm_elec = 0.544695364
        qq_elec = -0.009873526

    elif N_bp == 3:
        x_bp_val = np.array([0.069264723, 0.56492096, 1.09206815])
        y_bp_val = np.array([0.027854647, 0.312682424, 0.584970932])
        mm_elec = np.array([0.574647821, 0.51653222])
        qq_elec = np.array([-0.011948175, 0.020882546])

    elif N_bp == 4:
        x_bp_val = np.array(
            [0.069264723, 0.378789088, 0.735922253, 1.09206815])
        y_bp_val = np.array(
            [0.027854647, 0.210531814, 0.403454818,	0.584970932])
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
        

    return x_bp_val, y_bp_val, mm_elec, qq_elec

# ------------------------------------------------------------------------------
# FUEL CELL
# ------------------------------------------------------------------------------
"""
Reference: Morandi, Lisa (2022): Identification and MPC control of a hydrogen 
energy storage system for peak shaving EV demand. ETH Zurich. Online verfügbar 
unter https://www.research-collection.ethz.ch/handle/20.500.11850/590237.

FC Model: Input = Vdot_FC_H2 Flow, Output = P_FC

Vdot_FC_H2 = c * i_FC 
where    c = #cells * M / (zF * rho) * 3600s/h * 1000 Nl/Nm3 * 1 / (60min/h) * i_FC
           = #cells * 0.007 Nl/(A*min) * i_FC
        
For 100kW FC consisting of 2 stacks with a capacity of 50kW each:

#cells = 2 stacks * 238 cells/stack = 476 => Number of cells
"""

def pwa_FC():
    F        = 9.648533212331*10**4      # Faraday constant in [A*s/mol]
    M        = 22.414                    # molar volume of an ideal gas in [Nl/mol]
    z        = 2                         # Number of electrons exchanged in reaction

    N_stacks          = 2
    N_cells_per_stack = 238
    N_cells           = N_stacks * N_cells_per_stack

    c = N_cells *  M / (z * F) * (3600 / 1000)          # in [Nm3/Ah]
    
    FC_param_struct = {
        'P_FC_max'       : 100*1000,                         # Maximum power in [kW]
        'P_FC_min'       : 10*1000,                          # Minimum power in [kW]
        'FC_i_breakpoint': 123.054000168826,                 # Breakpoint current in A
        'FC_p_breakpoint': 47.96736518*1000,                 # Breakpoint power in Watts (converted from kW)
        'FC_slope_high_current': 3.31471893412921/1000,           # Slope above breakpoint in [1/kV]
        'FC_slope_low_current' : 2.56536917754678/1000,           # Slope below breakpoint in [1/kV]
        'FC_i_to_VH2' : c,
        'FC_dc_to_ac': 0.956183048244058                     # relation between dc input and ac power delivered
    }
    
    return FC_param_struct


# FC Model Validation

W_H2_HHV = 3.545*1000                       # The vol specific chemical energy (at STP) for HHV in [Wh/m3]

FC_param_struct = pwa_FC()

P_FC_out   = np.linspace(FC_param_struct['P_FC_min'], FC_param_struct['P_FC_max'])
i_FC       = [0] * len(P_FC_out)
Vdot_FC_H2 = [0] * len(P_FC_out)
P_FC_in    = [0] * len(P_FC_out)

for i in range(len(P_FC_out)):
    if FC_param_struct['P_FC_min'] <= P_FC_out[i] <= FC_param_struct['FC_p_breakpoint']:
        i_FC[i]       = FC_param_struct['FC_slope_low_current'] * P_FC_out[i]
        Vdot_FC_H2[i] = FC_param_struct['FC_i_to_VH2'] * i_FC[i]
        P_FC_in[i]    = W_H2_HHV * Vdot_FC_H2[i]
         
    elif FC_param_struct['FC_p_breakpoint'] < P_FC_out[i] <= FC_param_struct['P_FC_max']:
        i_FC[i]       = FC_param_struct['FC_slope_high_current'] * (P_FC_out[i] - FC_param_struct['FC_p_breakpoint']) + FC_param_struct['FC_i_breakpoint']
        Vdot_FC_H2[i] = FC_param_struct['FC_i_to_VH2'] * i_FC[i]
        P_FC_in[i]    = W_H2_HHV * Vdot_FC_H2[i]  # Corrected to index Vdot_FC_H2[i]
        
    else:
        P_FC_in[i] = 0
    

# Convert P_FC_out and P_FC_in to kW by dividing by 1000

P_FC_out_kW = np.array(P_FC_out) / 1000
P_FC_in_kW = np.array(P_FC_in) / 1000

# Calculate efficiency in %
efficiency = (P_FC_out_kW / P_FC_in_kW) * 100

# Create the figure and axes
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

# Plot i_FC vs Vdot_FC_H2 without colorbar
sc1 = ax1.scatter(Vdot_FC_H2, i_FC, c='blue')  # Using a solid color for simplicity
ax1.set_xlabel('Vdot_FC_H2', fontsize=20)
ax1.set_ylabel('i_FC', fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=20)

# Plot P_FC_in_kW vs P_FC_out_kW with colorbar
sc2 = ax2.scatter(P_FC_in_kW, P_FC_out_kW, c=efficiency, cmap='viridis')
ax2.set_xlabel('P_FC_in (kW)', fontsize=20)
ax2.set_ylabel('P_FC_out (kW)', fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
cb2 = fig.colorbar(sc2, ax=ax2)
cb2.set_label('Efficiency (%)', fontsize=20)
cb2.ax.tick_params(labelsize=20)
# ax2.set_ylim([0, 200])

# Tight layout to make sure labels don't overlap and everything fits well
plt.tight_layout()
plt.show()


# ------------------------------------------------------------------------------
# COMPRESSOR
# ------------------------------------------------------------------------------
"""
Sources: thesis (Josien de Koning and) Timo Laaksonlaita
"""


def compressor_curve():
    mm_c = np.array([4.8756, 6.5348, 2.5470])
    med_c = 4.078  # Midway H2 flow compressor [kg/h]
    return mm_c, med_c
