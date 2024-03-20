# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:55:43 2024

@author: fism
Battery model from Hanmin:

Experimental implementation of an emission-aware prosumer with online flexibility quantification and provision

URL: https://arxiv.org/abs/2110.12831

"""
import gurobipy as gp
from gurobipy import Model, GRB

#------------------------------------------------------------------------------
# Constants and parameters
#------------------------------------------------------------------------------

A_ebat     = 0.998           # Efficiency loss coefficient
B_ebat     = [0.256, 0.251]  # Matrix for power flow in/out of the battery
SOC_min    = 20              # Minimum state of charge for better lifetime
SOC_max    = 80              # Maximum state of charge for better lifetime
P_ebat_max = ...             # Maximum power for charging/discharging, to be defined by the user
nHours     = 24              # Number of hours, to be defined by the user

#------------------------------------------------------------------------------
# Model setup
#------------------------------------------------------------------------------

model = Model("Battery_Model")

#------------------------------------------------------------------------------
# Decision variables
#------------------------------------------------------------------------------
E_b          = model.addVars(nHours, name="E_b")                               # State of Charge
P_ch         = model.addVars(nHours, name="P_ch")                              # Power for charging
P_disch      = model.addVars(nHours, name="P_disch")                           # Power for discharging
epsilon_ebat = model.addVars(nHours, name="epsilon_ebat")                      # Slack variables
z_ebat       = model.addVars(nHours, vtype=GRB.BINARY, name="z_ebat")          # Binary variable for mutual exclusiveness

#------------------------------------------------------------------------------
# Constraints
#------------------------------------------------------------------------------

for t in range(nHours):
    # E_b update equation constraint (equation 11)
    if t < nHours - 1:  # No E_b[t+1] for the last hour
        model.addConstr(E_b[t+1] == A_ebat * E_b[t] + B_ebat[0] * P_disch[t] + B_ebat[1] * P_ch[t])

    # Charging power constraints (equation 12)
    model.addConstr(0 <= P_ch[t])
    model.addConstr(P_ch[t] <= P_ebat_max * z_ebat[t])

    # Discharging power constraints (equation 13)
    model.addConstr(-P_ebat_max * (1 - z_ebat[t]) <= P_disch[t])
    model.addConstr(P_disch[t] <= 0)

    # State of charge limits (equation 14)
    model.addConstr(SOC_min - epsilon_ebat[t] <= E_b[t])
    model.addConstr(E_b[t] <= SOC_max)

    # Slack variable constraint (equation 15)
    model.addConstr(epsilon_ebat[t] >= 0)

#------------------------------------------------------------------------------
# Objective function placeholder 
#------------------------------------------------------------------------------

# model.setObjective(...)

#------------------------------------------------------------------------------
# Optimization
#------------------------------------------------------------------------------

#model.optimize()

