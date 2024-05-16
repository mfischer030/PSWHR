%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                                 PEMEC                                   %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic model of a proton exchange membrane electrolyzer (PEMEC).
% A lumped, first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM electrolyzer are captured.

% The data represent a big scale PEM electrolyzer, up to 100 kW.
% It is modeled based on the commercial product Silyzer100 of Siemens
% (100 kW).
% The operative power can range in 5-100% of the rated power.
% The hydrogen and oxygen are assumed to be produced at 40 bar.
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
eta_p   = 0.7;           % Pump efficiency
tau_p   = 0.5;           % Pump time constant                              [s]

% Pressure variation from water pump [Pa]
DeltaP = (p_an - p_atm) / eta_p;

% RECTIFIER EFFICIENCY
eta_ACDC = 0.96;

%%                            PEMEC STACK 

% UTILIZATION FACTOR
c1   = [I_des * 0.05, I_des * 0.5, I_des]';
c2   = [0.99, 0.95, 0.90]';
polU = polyfit(c1, c2, 2);
U    = polyval(polU, I(1));

% MODELING PARAMETERS - Cell electrochemistry
CH        = 1000;       % Concentration of ion H+ (1200 or 1000)           [mol / m3]
cH2_me0   = 10.0;       % H2 reference concentration                       [mol / m3]
cO2_me0   = 0.10;       % O2 reference concentration                       [mol / m3]
DH        = 1.5e-09;    % Diffusion coefficient of ion H+                  [cm2 / s]
Epro      = 20000;      % Diffusion activation energy                      [J / mol]
io_an     = 3e-04;      % Exchange current density - anode                 [A / m2]
io_cat    = 1e+02;      % Exchange current density - cathode               [A / m2]
Req_el    = 180e-06;    % Equivalent electric resistance                   [ohm]
t_an      = 0.00020;    % Anode thickness                                  [m]
t_cat     = 0.00020;    % Cathode thickness                                [m]
t_mem     = 0.00020;    % Membrane thickness                               [m] 
epsilon   = 0.3;        % Electrode porosity
epsilon_p = 0.11;       % Electrode percolation threshold
lambda    = 25;         % Degree of humidification of the membrane         [molH2O / molSO3-]

% FARADAY'S EFFICIENCY DATA - Data from [1,4]
f1 = 25000;             % Faraday coefficient                              [A2 / m4]
f2 = 0.98;              % Faraday coefficient      

% THERMAL PARAMETER
C_th   = 250e03;        % Thermal capacitance from [2]                     [J / K]
h_cond = 380;           % Conduction coefficient                           [W / C]
h_conv = 0.1;           % Convection coefficient                           [W / C A] 
R_th   = 0.0817;        % Thermal resistance from [2]                      [K / W]

% Cooling system
c_cw   = 4.18e03;       % Cooling water specific heat                      [J / kg C]
m_cw   = 0.250000;      % Cooling water mass flow                          [kg / s]
C_cw   = c_cw * m_cw;   % Cooling water thermal capacity                   [W / C]
T_cwin = 20 + 0;        % Inlet cooling water temperature                  [C]
V_tn   = 1.48;          % Thermo-neutral voltage                           [V]
tau_th = C_th * R_th;   % Thermal time constant                            [s]

% INITIAL PARTIAL PRESSURES [atm]
p_H20  = y_H2O * p_an / 101325;
p_H2O0 = y_H2 * p_an / 101325;
p_O20  = y_O2 * p_cat / 101325;

% DIFFUSION OVERPOTENTIAL COEFFICIENTS [V]
a = 3.64e-04;   % Coefficient from Marangio et al. (2009)
b = 2.334;      % Coefficient from Marangio et al. (2009)

% H2-H2O diffusion coefficient [cm2 / s]
D_H20 = ((33.3 * 647.3)^(1 / 3)) * ((12.8 * 218.3) ^ (5/12)) * ...
    ((1/2 + 1/18) ^ (1/2)) / (p_cat*1e-05);  

% O2-H2O diffusion coefficient [cm2 / s]
D_O20 = ((154.4 * 647.3) ^ (1/3)) * ((49.7 * 218.3) ^ (5/12)) * ...
    ((1/32 + 1/18) ^ (1/2)) / (p_an*1e-05);  

% H2-H2O effective diffusion coefficient  [cm2 / s]
Deff_cat_ = D_H20 * epsilon * (((epsilon - epsilon_p) / ...
    (1 - epsilon_p)) ^ 0.785); 

% O2-H2O effective diffusion coefficient  [cm2 / s]
Deff_an_  = D_O20 * epsilon * (((epsilon - epsilon_p) / ...
(1 - epsilon_p)) ^ 0.785);   
  
% Unit conversion - cm2 ---> m2
Deff_an0  = Deff_an_  * 1e04;   
Deff_cat0 = Deff_cat_ * 1e04; 

% ACTIVATION OVERPOTENTIAL [V]
alpha_an  = 2;                      % Anode charge tranfer coefficient
alpha_cat = 0.5;                    % Cathode charge tranfer coefficient

% Calculation parameters
K2 = 1 / (2 * A_cell * io_an);
K4 = 1 / (2 * A_cell * io_cat);

%%                     PEMEC - H2 PRODUCTION MODEL 

% WATER TRANSPORT
% Constant parameters - Values from [3]
Dw = 1.28e-10;              % Water diffusion coefficient                  [m2 / s]
Kw = 1.58e-18;              % Darcy's constant                             [m2]
nd = 7;                     % Electro-osmotic drag coeff.                  [molH2O / molH+]

% Water transport for pressure gradient [mol / s]
n_H2Opg = Kw * A_cell * rho_H2O / (mu_H2O * M_H2O * t_mem)*(p_cat - p_an);

% Overall water through the membrane [mol / s]
n_H2Omem2 = 1 - Dw / t_mem * (t_cat / Deff_cat_ + t_an / Deff_an_);
n_H2Omem1 = Dw / t_mem * (t_an/(2 * F * Deff_an_));

%%                          DYNAMIC BEHAVIOR

% ANODE AND CATODE CHANNELS - Data from [3]
% Bipolar capacity [F]
C  = 0.008; 

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
p_H2Oout = y_H2Oout * p_an / 101325;       % H2O partial pressure          [atm]
p_O2out  = y_O2out * p_an / 101325;        % O2 partial pressure           [atm]
p_H2out  = y_H2out * p_cat / 101325;       % H2 partial pressure           [atm]

% Partial pressure at inlet
p_H2Oin = p_H2O0;                          % H2O partial pressure          [atm]
p_O2in  = p_H2O0;                          % O2 partial pressure           [atm]
p_H2in  = p_H2O0;                          % H2 partial pressure           [atm]

% Valve constants
K_H2_  = n_H2r / p_H2out / 101325;         % H2 valve constant             [mol / Pa s]   
K_H2O_ = n_H2Oout / p_H2Oout / 101325;     % H2O valve constant            [mol / Pa s]
K_O2_  = n_O2r / p_O2out / 101325;         % O2 valve constant             [mol / Pa s] 
K_H2   = n_H2r / p_H2out;                  % H2 valve constant             [mol / atm s]   
K_H2O  = n_H2Oout / p_H2Oout;              % H2O valve constant            [mol / atm s]
K_O2   = n_O2r / p_O2out;                  % O2 valve constant             [mol / atm s] 

% Time constants
tau_H2  = 2.07;                            % H2 time constant              [s]
tau_H2O = 4.41;                            % H2O time constant             [s]
tau_O2  = 4.14;                            % O2 time constant              [s]

% Calculation parameters [mol/As]
K_r = N_cell / (4 * F);         

%%                              STORAGE

% DATA FROM PSI PLANT  
VS_H2 = 40;                             % Volume of hydrogen storage       [m3]
VS_O2 = 20;                             % Volume of oxygen storage         [m3]
p_max = 45e05;                          % Maximum pressure of the storage  [Pa]
zeta  = 1;                              % Hydrogen compressibility factor
