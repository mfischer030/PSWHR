%--------------------------------------------------------------------------------------------------%
%                                      H2-O2 PEM Fuel Cell                                         %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                           Process Engineering Institute, June 2016


%--------------------------------------------------------------------------------------------------%
%                                   Stack and BoP Simulation Data                                  %
%--------------------------------------------------------------------------------------------------%

%%                               DESCRIPTION AND ASSUMPTIONS

% Thermodynamic model of a H2-O2 proton exchange membrane fuel cell (PEMFC). 
% A lumped first-principle approach is followed. Static and dynamic behavior are described, and both 
% electrochemical and thermal features of the PEM fuel cell are captured.

% The model is based on the system Belenos/PSI B25 (25 kw) using H2 as a fuel and pure oxygen as an
% oxidant in order to increase stack performance and power density [8]. 
% A first order dynamics is used to simulate the time response of both the electrical and thermal 
% produced powers. The time constant for the electrical power is typically 1 order of magnitude 
% smaller than that of the thermal power.

%%                                 DESIGN SPECIFICATIONS

% TEMPERATURES
T_amb = T_a + 273.15;     % Ambient temperature                            [K]
T_des = T - 273.15;       % Design temperature                             [C]
T0    = T_a;              % Cell initial temperature                       [C]

% PRESSURES OF STORAGE
F_pst = 0;                % Storage pressure factor

%%                                   THERMODYNAMIC DATA

% MOLAR WEIGHT
M_H2  = 2;                % H2 molar weight                   [kg/kmol]
M_H2O = 18;               % H2O molar weight                  [kg/kmol]
M_O2  = 32;               % N2 molar weight                   [kg/kmol]

% LOWER HEATING VALUE
LHV_H2 = 119.96e03;       % H2 lower heating value            [J/g]

% JOULE-THOMPSON COEFFICIENT
JT_H2 = -0.038;           % H2 JT coefficient (20 bar, 70 °C) [K/bar]
JT_O2 = 0.194;            % O2 JT coefficient (20 bar, 70 °C) [K/bar]

% SPECIFIC VOLUME
mu_H2  = 0.0116;          % H2 specific volume                [m3/mol]                  
mu_O2  = 0.0115;          % O2 specific volume                [m3/mol]
mu_H2O = 1.8468e-05;      % H2O specific volume               [m3/mol]

%%                                    BALANCE OF PLANT

% DC/AC CONVERTER
eta_DCAC = 0.95;          % conversion efficiency of AC/DC converter
eta_ph   = 0.97;          % Pre-heater efficiency

% GAS PRE-HEATER
% H2 specific pre-heating enthalpy
hH2_ph = (refpropm('H', 'T', T, 'P', p_an * 100, 'HYDROGEN.FLD') - ...
  refpropm('H', 'T', T_amb, 'P', p_an * 100, 'HYDROGEN.FLD')) * M_H2 / 1000;

% H2 specific heat [J/molK]
CH2 = refpropm('C', 'T', T_amb, 'P', p_an * 100, 'HYDROGEN.FLD') * M_H2 / 1000;

% O2 specific heat [J/molK]
CO2 = refpropm('C', 'T', T_amb, 'P', p_cat * 100, 'OXYGEN.FLD') * M_O2 / 1000;

% O2 specific pre-heating enthalpy [J/mol]
hO2_ph = (refpropm('H', 'T', T, 'P', p_cat * 100, 'OXYGEN.FLD') - ...
  refpropm('H', 'T', T_amb, 'P', p_cat * 100, 'OXYGEN.FLD')) * M_O2 / 1000;

% MISCELLANEOUS AUXILIARIES
P_gamma = 100.0;            % BoP fixed consumption                   [W]

% GAS RELATIVE HUMIDITY
% Water saturation pressure [atm]
pH2Osat = WaterSaturationPressure (T);

% Water molar fraction in fed gas
xH2O_H2 = rH_H2 * pH2Osat / p_an;
xH2O_O2 = rH_O2 * pH2Osat / p_cat;

% H2O specific pre-heating enthalpy
hH2O_ph = (refpropm('H', 'T', T, 'P', p_an * 100, 'WATER.FLD') - ...
  refpropm('H', 'T', T_amb, 'P', p_an * 100, 'WATER.FLD')) * M_H2O / 1000;

%%                                  PEMFC STACK

% MASS BALANCES 
% Inlet molar flow [mol/s]
nH2in     = I * A_cell * N_cell / (2 * F) * lambda_H2;
nO2in     = I * A_cell * N_cell / (4 * F) * lambda_O2;
nH2Oin_H2 = nH2in * xH2O_H2 / (1 - xH2O_H2);
nH2Oin_O2 = nO2in * xH2O_O2 / (1 - xH2O_O2);

% Reacting molar flow [mol/s]
nH2r     = I * A_cell * N_cell / (2 * F);
nO2r     = I * A_cell * N_cell / (4 * F);
nH2Or_H2 = I * A_cell * N_cell / (2 * F);
nH2Or_O2 = 0;

% Outlet molar flow [mol/s]
nH2out     = nH2in - nH2r;
nO2out     = nO2in - nO2r;
nH2Oout_H2 = nH2Oin_H2 + nH2Or_H2;
nH2Oout_O2 = nH2Oin_O2;

% Inlet partial pressure [atm]
p_H20  = (1 - xH2O_H2) * p_an;
p_O20  = (1 - xH2O_O2) * p_cat;
p_H2O0 = xH2O_H2 * p_an;

% Outlet partial pressure [atm]
pH2out     = nH2out / (nH2out + nH2Oout_H2) * p_an;
pO2out     = nO2out / (nO2out + nH2Oout_O2) * p_cat;
pH2Oout_H2 = nH2Oout_H2 / (nH2out + nH2Oout_H2) * p_an;

% Average partial pressure [atm]
pH2avg     = (p_H20 + pH2out) * 0.35;
pO2avg     = (p_O20 + pO2out) * 0.35;

% ELECTROCHEMICAL MODEL
% Valve constants
K_H2   = pH2out / nH2out;           % H2 valve constant  [mol/atm s] 
K_O2   = pO2out / nO2out;           % O2 valve constant  [mol/atm s] 
K_H2O  = pH2Oout_H2 / nH2Oout_H2;   % H2O valve constant [mol/atm s]

% Time constants
tau_H2  = 0.27;     % H2 time constant   [s]
tau_H2O = 0.81;     % H2O time constant  [s]
tau_O2  = 0.47;     % O2 time constant   [s]

% Calculation parameters
K_r = N_cell / (4 * F);   % Reaction parameter [mol/C]  

% Diffusion losses - Data from [2]
B     = 0.005;    % fitting coefficient for concentration losses [V]      
i_max = 2.027;    % maximum current density                      [A/cm2] 

% Activation losses - Data from [4]
c_O2 = pO2avg / (5.08e06 * exp(-498 / T));   % O2 concentration at TPB [mol / cm3]

% Fitting parameters - Data from [1]
xi1 = -0.61693;
xi2 = 0.0016838;
xi3 = 1.9385e-06;
xi4 = -6.4922e-05;

% Ohmic losses - Data from [4]
R_el = 0.00025492;   % equivalent internal resistance [ohm] 
psi  = 23;           % membrane wettability

% THERMAL MODEL
C_th   = 350e03;        % Thermal capacitance from [2]    [J / K]
c_cw   = 4.18e03;       % Cooling water specific heat     [J / kg C]
m_cw   = 0.25000;       % Cooling water mass flow         [kg / s]
C_cw   = c_cw * m_cw;   % Cooling water thermal capacity  [W / C]
h_cond = 500;           % Conduction coefficient          [W / C]
h_conv = 0.1;           % Convection coefficient          [W / C A] 
R_th   = 0.0817;        % Thermal resistance from [2]     [K / W]
T_cwin = 20;            % Inlet cooling water temperature [C]
V_tn   = 1.253;         % Thermo-neutral voltage          [V]
tau_th = C_th * R_th;   % Thermal time constant           [s]

%%                            LIBRARIES

% WATER SATURATION PRESSURE
chi1 = 0.00129697;   % Fitting coefficient - Thermodynamic
chi2 = - 1.529053;   % Fitting coefficient - Thermodynamic
chi3 = 681.731481;   % Fitting coefficient - Thermodynamic
chi4 = - 136025.0;   % Fitting coefficient - Thermodynamic
chi5 = 10234070.7;   % Fitting coefficient - Thermodynamic
