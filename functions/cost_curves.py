# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:50:58 2024

@author: fism
"""

import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Heat Exchanger Cost Curve
#------------------------------------------------------------------------------

A = np.linspace(0.5, 30, 100)
c = np.log(A) - 0.6395 * A**2 + 947.2 * A + 227.9

N_bp   = 3
idx    = np.round(np.linspace(0, 99, N_bp)).astype(int)
x_val  = A[idx]
y_val  = c[idx]
a_etae = (y_val[1:] - y_val[:-1]) / (x_val[1:] - x_val[:-1])
b_etae = y_val[1:] - a_etae * x_val[1:]
aa_HEX = a_etae
bb_HEX = b_etae

#------------------------------------------------------------------------------
# Compressor Cost 
#------------------------------------------------------------------------------

S_c = np.linspace(5, 20, 100)
C_c = 43872 * S_c**0.5861

xx_c = np.array([5, 12.5, 20])
yy_c = 43872 * xx_c**0.5861
aa_c = (yy_c[1:] - yy_c[:-1]) / (xx_c[1:] - xx_c[:-1])
bb_c = yy_c[:-1] - aa_c * xx_c[:-1]

#------------------------------------------------------------------------------
# Electrolyser Cost Curve 
#------------------------------------------------------------------------------
"""A. H. Reksten, M. S. Thomassen, S. M˜ A¸ller-Holst, and K. Sundseth, 
“Projecting the future cost of pem and alkaline water electrolysers; 
a capex model including electrolyser plant size and technology development,” 
International Journal of Hydrogen Energy, vol. 47, pp. 38106–38113, 11 2022."""

V_0   = 2020     # reference year => from paper
V     = 2023     # plant installation year            => from paper
alpha = 0.622    # fitting constant / scaling factor  => from paper
beta  = -158.9   # fitting constant / learning factor => from paper
k_0   = 585.85   # fitting constant                   => from paper
k     = 9458.2   # fitting constant                   => from paper
SE    = 510      # Standard error                     => from paper

S_ELY = np.linspace(0.01, 2000, 2000)                                # Compressor size (kW)
C_e = (k_0 + k / S_ELY * S_ELY**alpha) * (V / V_0)**beta * S_ELY     # Compressor cost

xx_e = np.array([250, 400, 550])
yy_e = (k_0 + k / xx_e * xx_e**alpha) * (V / V_0)**beta * xx_e
aa_e = (yy_e[1:] - yy_e[:-1]) / (xx_e[1:] - xx_e[:-1])
bb_e = yy_e[:-1] - aa_e * xx_e[:-1]

#------------------------------------------------------------------------------
# Plotting Cost Curves
#------------------------------------------------------------------------------

# Updating font sizes for legends, axis values, and labels to 24pt

fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Customizing fonts
legend_fontsize = 24
axis_label_fontsize = 24
axis_values_fontsize = 24
title_fontsize = 24

# Heat Exchanger Cost Curve with cost in kEUR
axs[0].plot(A, c / 1000, label='Cost curve')  # Dividing cost by 1000
for i in range(len(aa_HEX)):
    x_range = A[idx[i]:idx[i+1]+1] if i < len(aa_HEX) - 1 else A[idx[i]:]
    axs[0].plot(x_range, (aa_HEX[i]*x_range + bb_HEX[i]) / 1000, label=f'PWA pt{i+1}')  # Dividing cost by 1000
axs[0].legend(fontsize=legend_fontsize)
axs[0].set_xlabel('Heat exchanger area [m^2]', fontsize=axis_label_fontsize)
axs[0].set_ylabel('Cost [kEUR]', fontsize=axis_label_fontsize)  # Adjusted to kEUR
axs[0].set_title('Heat Exchanger Cost Curve', fontsize=title_fontsize)
axs[0].tick_params(axis='both', which='major', labelsize=axis_values_fontsize)

# Compressor Cost Curve with cost in kEUR
axs[1].plot(S_c, C_c / 1000, label='Actual curve')  # Dividing cost by 1000
axs[1].plot(S_c[:67], (aa_c[0]*S_c[:67] + bb_c[0]) / 1000, label='PWA pt1')  # Dividing cost by 1000
axs[1].plot(S_c[66:], (aa_c[1]*S_c[66:] + bb_c[1]) / 1000, label='PWA pt2')  # Dividing cost by 1000
axs[1].legend(fontsize=legend_fontsize)
axs[1].set_xlabel('Compressor size [kW]', fontsize=axis_label_fontsize)
axs[1].set_ylabel('Compressor cost [kEUR]', fontsize=axis_label_fontsize)  # Adjusted to kEUR
axs[1].set_title('Compressor Cost Curve', fontsize=title_fontsize)
axs[1].tick_params(axis='both', which='major', labelsize=axis_values_fontsize)

# Electrolyser Cost Curve
axs[2].plot(S_ELY, C_e / 1000, label='Cost curve')  # Already in kEUR
axs[2].plot(S_ELY[280:401], (aa_e[0]*S_ELY[280:401] + bb_e[0])/1000, label='PWA pt1')  # No change needed
axs[2].plot(S_ELY[400:511], (aa_e[1]*S_ELY[400:511] + bb_e[1])/1000, label='PWA pt2')  # No change needed
axs[2].set_xlim([100, 150])
axs[2].set_ylim([200, 250])
axs[2].legend(fontsize=legend_fontsize)
axs[2].set_xlabel('Electrolyser size [kW]', fontsize=axis_label_fontsize)
axs[2].set_ylabel('Electrolyser cost [kEUR]', fontsize=axis_label_fontsize)  # Already correctly labeled
axs[2].set_title('Electrolyser Cost Curve', fontsize=title_fontsize)
axs[2].tick_params(axis='both', which='major', labelsize=axis_values_fontsize)

plt.tight_layout()
plt.show()
