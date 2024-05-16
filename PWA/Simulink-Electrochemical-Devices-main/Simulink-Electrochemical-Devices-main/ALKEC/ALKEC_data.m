%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                                 ALKEC                                   %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic model of an alkaline electrolyzer (ALKEC).
% A lumped, first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the alkaline electrolyzer are captured.

% The data represent a big scale alkaline electrolyzer, up to 100 kW.
% It is modeled based on the commercial product HySTAT10 of Hydrogenics
% (54 kW).
% The operative power can range in 15-100% of the rated power.
% The hydrogen is assumed to be produced at 10 bar and stored at 
% (p, T) = (10 bar, 25 °C).
% A reference temperature of 70 °C is chosen as the commercial product is
% assumed to be controlled to work close to the design temperature in the
% whole range of power.
% For the moment, heat integration is not considered for the electrolyzer,
% although the device can produce heat while producing hydrogen at high
% power.
% A first order dynamics is used to simulate the time response of both the
% electrical and thermal produced powers. The time constant for the
% electrical power is typically 1 order of magnitude smaller than that of
% the thermal power.
 
% Modeling parameter were chosen from the literature. The following
% assumptions were made:
% 1. 0-D model: greybox approach;
% 2. Ideal gases;
% 3. Negligible pressure drop along the channels;
% 4. Constant operating temperature (design conditions).

%%                        DESIGN SPECIFICATIONS

% TEMPERATURES
T_amb = T_a + 273.15;     % Ambient temperature                            [K]
T_an  = T;                % Anode temperature                              [K]               
T_cat = T;                % Cathode temperature                            [K] 
T_des = T - 273.15;       % Design temperature                             [C]
T_max = 80;               % Maximum temperature of the EC                  [C]
T0    = T_a;              % Cell initial temperature                       [C]          

%%                            THERMODYNAMICS

% MOLAR WEIGHT [kg/kmol]
M_H2    = 2.01588;          % Molar mass of H2                             [g / mol] 
M_H2O   = 18.01528;         % Molar mass of H2O                            [g / mol]
M_O2    = 31.9988;          % Molar mass of O2                             [g / mol]

% MASS DENSITY
v_H2    = 0.08988e03;       % Volumic mass of H2 - gas                     [g / m^3]
v_H240  = 3.1777e03;        % Volumic mass of H2 at (p,T) = (40 bar, 25 C) [g / m^3]
v_H250  = 3.4443e03;        % Volumic mass of H2 at (p,T) = (50 bar, 70 C) [g / m^3]
v_H2O   = 999.7e03;         % Volumic mass of H2O - liquid                 [g / m^3]
v_O2    = 1.308e03;         % Volumic mass of O2 - gas                     [g / m^3]
v_O240  = 52.876e03;        % Volumic mass of H2 at (p,T) = (40 bar, 25 C) [g / m^3]
v_O250  = 56.761e03;        % Volumic mass of H2 at (p,T) = (50 bar, 70 C) [g / m^3]
rho_H2O = v_H2O/1000;       % Water density                                [g / m3]
mu_H2O  = 1.1e-03;          % Water dynamic viscosity                      [Pa s]

%%                             AUXILIARES

% WATER PUMP
C_H2O   = 4.182;         % H2O specific heat at (p, T) = (1 bar, 323 K)    [kJ / kg K]
P_gamma = 2.5e03;        % Additional auxiliary power (gas cleaning)       [W]
eta_p   = 0.4;           % Pump efficiency
tau_p   = 0.5;           % Pump time constant                              [s]

% Pressure variation from water pump [Pa]
DeltaP = (p_an - p_atm) / eta_p;

% RECTIFIER EFFICIENCY
eta_ACDC = 0.95;

%%                            PEMEC STACK 

% UTILIZATION FACTOR
c1   = [I_des * 0.15, I_des * 0.5, I_des]';
c2   = [0.98, 0.94, 0.90]';
polU = polyfit(c1, c2, 2);
U    = polyval(polU, I(1));

% MODELING PARAMETERS - Cell electrochemistry
t_an      = 0.00010;       % Anode thickness                               [m]
t_cat     = 0.00010;       % Cathode thickness                             [m]
t_mem     = 0.00010;       % Membrane thickness                            [m] 

% MODELING PARAMETERS - Experimental Data [1]
% Overvoltage parameters
r1    = 3.53855e-04;       % Resistivity coefficient                       [ohm m2]
r2    = -3.02150e-06;      % Resistivity coefficient                       [ohm / m2 C]
s     = 0.22396;           % Activation overvoltage coefficient            [V]
t1    = 2.73093;           % Diffusion overvoltage coefficient             [m2 / A]
t2    = -2.40447e02;       % Diffusion overvoltage coefficient             [m2 C / A]
t3    = 5.99576e03;        % Diffusion overvoltage coefficient             [m2 C2 / A]
theta = 0.95;

% Faraday efficiency parameters - Data from [1]
f1 = 2.50;                 % Faraday coefficient                           [A2 / m4]
f2 = 0.97;                 % Faraday coefficient      

% THERMAL PARAMETER
C_th   = 625e03;           % Thermal capacitance                           [J / C]
h_cond = 250;              % Conduction coefficient                        [W / C]
h_conv = 0.97;             % Convection coefficient                        [W / C A] 
R_th   = 0.167;            % Thermal resistance                            [K / C]

% Cooling system
c_cw   = 4.18e03;          % Cooling water specific heat                   [J / kg C]
m_cw   = 0.250000;         % Cooling water mass flow                       [kg / s]
C_cw   = c_cw * m_cw;      % Cooling water thermal capacity                [W / C]
T_cwin = 20 + 0;           % Inlet cooling water temperature               [C]
V_tn   = 1.48;             % Thermo-neutral voltage                        [V]
tau_th = C_th * R_th;      % Thermal time constant                         [s]

% PARTIAL PRESSURES CALCULATION [Pa]
% Initial partial pressure [Pa]
p_H20  = y_H2O * p_an;
p_H2O0 = y_H2 * p_an;
p_O20  = y_O2 * p_cat;

% Partial pressure at inlet
p_H2Oin = p_H2O0;                          % H2O partial pressure          [atm]
p_O2in  = p_O20;                           % O2 partial pressure           [atm]
p_H2in  = p_H20;                           % H2 partial pressure           [atm]

% Molar flow rates 
n_H2r    = N_cell * I(1) / (2 * F);        % H2 mole flow reacting         [mol / s]
n_H2Or   = N_cell * I(1) / (2 * F);        % H2O mole flow reacting        [mol / s]
n_H2Oin  = N_cell * I(1) / (2 * F * U);    % H2O mole flow in              [mol / s]
n_H2Oout = n_H2Oin - n_H2Or;               % H2O mole flow out             [mol / s]
n_O2r    = N_cell * I(1) / (4 * F);        % O2 mole flow reacting         [mol / s]

% Total mole flow at outlet
n_anout = n_H2Oout + n_O2r;                % Anode mole flow               [mol / s]
n_ctout = n_H2r;                           % Cathode mole flow             [mol / s]

% Mole fractions at outlet
y_H2Oout = n_H2Oout / n_anout;             % H2O mole fraction - anode
y_O2out  = 1 - y_H2Oout;                   % O2 mole fraction - anode
y_H2out  = 1;                              % H2 mole fraction - cathode

% Partial pressure at outlet
p_H2Oout = y_H2Oout * p_an;                % H2O partial pressure          [atm]
p_O2out  = y_O2out * p_an;                 % O2 partial pressure           [atm]
p_H2out  = y_H2out * p_cat;                % H2 partial pressure           [atm]

% Average partial pressures
p_H2Oavg = (p_H2Oin + p_H2Oout) / 2;       % H2O partial pressure          [atm]
p_O2avg  = (p_O2in + p_O2out) / 2;         % O2 partial pressure           [atm]
p_H2avg  = (p_H2in + p_H2out) / 2;         % H2 partial pressure           [atm]

% Time constants
tau_E   = 2.00;                            % H2 time constant              [s]

%%                              STORAGE

% DATA FROM PSI PLANT  
VS_H2 = 40;                             % Volume of hydrogen storage       [m3]
VS_O2 = 20;                             % Volume of oxygen storage         [m3]
p_max = 45e05;                          % Maximum pressure of the storage  [Pa]
zeta  = 1;                              % Hydrogen compressibility factor
