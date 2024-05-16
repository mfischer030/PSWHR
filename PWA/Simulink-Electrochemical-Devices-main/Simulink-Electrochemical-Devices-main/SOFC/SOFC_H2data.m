%-------------------------------------------------------------------------%
%                     SOFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                             Process Engineering Institute, September 2015


%-------------------------------------------------------------------------%
%                          SOFC-CHP Parameters                            %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic modeling of a solid oxide fuel cell (SOFC). 
% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM fuel cell are captured.

% The model is based on the commercial product BlueGen of SolidPower, using
% methane as a fuel (1.6 kW). The polarization curve is also based on the
% same product.
% The operative power can range in 0-120% of the rated power.
% The fuel can use either natural gas or hydrogen as a fuel. The operation
% pressure would be higher when using hydrogen, assuming to have it in a 
% storage at a pressure higher than the ambient pressure.
% A reference temperature of 750 °C is chosen as the commercial product is
% assumed to be controlled to work close to the design temperature in the
% whole range of power.
% A first order dynamics is used to simulate the time response of both the
% electrical and thermal produced powers. The time constant for the
% electrical power is typically 1 order of magnitude smaller than that of
% the thermal power.
 
% Modeling parameter were chosen from the literature. The following
% assumptions were made:
% 1. 0-D model: greybox approach;
% 2. Ideal gases;
% 3. Dynamic describing a H2-fuelled channel;
% 4. Negligible pressure drop along the channels;
% 5. Constant operating temperature (design conditions).

%%                        DESIGN SPECIFICATIONS

% TEMPERATURES
T_amb = T_a + 273.15;     % Ambient temperature                            [K]
T_des = T - 273.15;       % Design temperature                             [C]
T0    = T_a;              % Cell initial temperature                       [C]

%%                       THERMODYNAMIC DATA

% MOLAR WEIGHT
M_air  = 28.97;           % air molar weight                               [kg / kmol]
M_Ar   = 40;              % Ar molar weight                                [kg / kmol]
M_CH4  = 16;              % CH4 molar weight                               [kg / kmol]
M_CO   = 28;              % CO molar weight                                [kg / kmol]
M_CO2  = 44;              % CO2 molar weight                               [kg / kmol]
M_C2H6 = 30.07;           % C2H6 molar weight                              [kg / kmol]
M_C3H8 = 44.10;           % C3H8 molar weight                              [kg / kmol]
M_H2   = 2;               % H2 molar weight                                [kg / kmol]
M_H2O  = 18;              % H2O molar weight                               [kg / kmol]
M_N2   = 28;              % N2 molar weight                                [kg / kmol]
M_O2   = 32;              % O2 molar weight                                [kg / kmol]         

% LOWER HEATING VALUE
LHV_CH4 = 50e03;          % CH4 lower heating value                        [kJ / kg]
LHV_CO  = 32e03;          % CO lower heating value                         [kJ / kg]
LHV_H2  = 120e03;         % H2 lower heating value                         [kJ / kg]

% SPECIFIC HEAT RATIO
gamma_air = 1.33;         % Air specific heat ratio
gamma_CH4 = 1.32;         % CH4 specific heat ratio
gamma_H2  = 1.409;        % H2 specific heat ratio

% SPECIFIC HEAT AT CONSTANT PRESSURE [J / g K]
% Prereformer conditions   
C_H2phout = 15;          % pre-reformer outlet - T = 1020 K   

%%                           AUXILIARIES

% EFFICIENCIES
eta_ACDC = 0.95;          % Electrical efficiency of AC/DC converter
eta_ph   = 0.92;          % Pre-heater efficiency

% TIME CONSTANTS
tau_c  = 0.10;            % Compressor time constant                       [s]
tau_ph = 5.00;            % Air pre-heater time constant                   [s]

% COMPRESSORS
beta_c = p_atm / p_an;    % Pressure ratio
eta_c  = 0.75;            % Compressor efficiency

% CH4 specific compression work [J / mol]  
h_H2c  = gamma_H2 * R * T / (gamma_H2 - 1) * ((1 / beta_c)^...
    ((gamma_H2 - 1) / gamma_H2) - 1) / eta_c;   

% Air specific compression work - air [J / mol]  
h_airc = gamma_air * R * T / (gamma_air - 1) * ((1 / beta_c)^...
    ((gamma_air - 1) / gamma_air) - 1) / eta_c; 

% MISCELLANEOUS AUXILIARIES
P_gamma = 100.0;          % Miscellaneous losses                           [W]
Q_gamma = 30.0;           % Heat-up losses                                 [W]

%%                           PRE-HEATING

% Inlet temperatures [K]
Tairin = 25 + 273.15;        % air inlet temperature                       [K]  
TH2in  = 15 + 273.15;        % H2O inlet temperature                       [K]          

% Outlet temperatures [K]
Tairout = T;                 % air outlet temperature                      [K]  
TH2out  = T;                 % fuel outlet temperature                     [K]  

% Inlet enthalpy [J/g]
hairph_in = refpropm('H', 'T', Tairin, 'P', p_cat/1000, 'AIR.PPF')/1000;
hH2ph_in  = refpropm('H', 'T', TH2in, 'P', p_an/1000, 'HYDROGEN.FLD')/1000;

% Outlet enthalpy [J/g]
hairph_out = refpropm('H', 'T', Tairout, 'P', p_cat/1000, 'AIR.PPF')/1000;

% Prereformer enthalpies [J/mol]
h_airph = (hairph_out - hairph_in) * M_air / eta_ph;
h_H2ph  = (C_H2phout * TH2out - hH2ph_in) * M_H2 / eta_ph;

%%                           AIR TREATMENT

% AIR SPLIT FACTOR (Combustor-to-FC)
F_airsplit = 1.1;

% AIR INLET FLOW - Aspen simulations
n_air     = n_H2in(1) * lambda;                % Air inlet to FC           [mol/s]
n_aircomb = n_H2in(1) * lambda * F_airsplit;   % Air inlet to combustor    [mol/s]

%%                            PEMFC STACK

% UTILIZATION FACTOR
c1   = [n_fueldes(1) * 0.1, n_fueldes(1), n_fueldes(1) * 1.25]';
c2   = [0.99, 0.90, 0.85]';
c3   = linspace(c1(1), c1(end), 50);
polU = polyfit(c1, c2, 2);
U    = polyval(polU, n_H2in(1));

% ELECTROCHEMICAL DATA
A_an      = 1.05e-05;       % Anode ohmic loss coefficient [1]             [ohm cm]
A_cat     = 2.38e-05;       % Cathode ohmic loss coefficient [1]           [ohm cm]
A_el      = 0.03 / T;       % Electrolite ohmic loss coefficient [1]       [ohm cm]
B_an      = -1150;          % Anode ohmic loss coefficient [1]             [K]
B_cat     = 1200;           % Cathode ohmic loss coefficient [1]           [K]
B_el      = 10300;          % Electrolite ohmic loss coefficient [1]       [K]
E_act_an  = 75e03;          % Anode activation energy [3]                  [J / mol]
E_act_cat = 96e03;          % Cathode activation energy [3]                [J / mol]
i_max     = 0.88;           % Limiting current density                     [A / cm2]
t_an      = 0.0545;         % Anode thickness [1]                          [cm]
t_cat     = 0.004;          % Cathode thickness [1]                        [cm]
t_el      = 0.0005;         % Electrolyte thickness [1]                    [cm]
alpha     = 0.50;           % Charge transfer coefficient [2]
gamma_an  = 5.5e4;          % Anode activation overpotential factor [3]    [A / cm2]
gamma_cat = 7.0e4;          % Cathode activation overpotential factor [3]  [A / cm2]
theta     = 0.910;          % OCV experimental deviation from theory [1]   

% MASS BALANCES
% Inlet molar fractions
y_H2in  = 0.995;
y_H2Oin = 1 - y_H2in;
x_O2in  = x_O2;

% Cell current [A]
n_fuel = n_H2in(1);
I = n_fuel * 2 * F * U / N_cell; 

% Calculation parameters
K_r = N_cell / (4 * F);              % Reaction parameter                  [mol / C]
Kr  = I * K_r;

% Inlet molar flows [mol/s]
n_H2in_  = n_fuel * y_H2in;
n_H2Oin = n_fuel * y_H2Oin;
n_O2in  = n_air * x_O2in;

% Outlet molar flows [mol/s]
n_H2out  = n_H2in_ - 2 * Kr;
n_H2Oout = n_H2Oin + 2 * Kr;
n_O2out  = n_O2in - Kr;

% Inlet partial pressures [atm]
p_H20  = y_H2in * p_an / 101325;
p_H2O0 = y_H2Oin * p_an / 101325;
p_O20  = x_O2in * p_cat / 101325;
                
% ANODE AND CATHODE CHANNELS DYNAMIC
% Time and valve constants
tau_H2  = 26.1;                      % H2 time constant                    [s]
tau_H2O = 78.30;                     % H2O time constant                   [s]
tau_O2  = 2.91;                      % O2 time constant                    [s]
K_H2    = 8.34e-01;                  % H2 valve constant                   [mol / atm s]   
K_H2O   = 2.81e-01;                  % H2O valve constant                  [mol / atm s]
K_O2    = 2.52e-00;                  % O2 valve constant                   [mol / atm s]

% THERMAL MODEL
C_th   = 500e03;            % Thermal capacitance from [2]                 [J / K]
c_cw   = 4.18e03;           % Cooling water specific heat                  [J / kg C]
m_cw   = 0.25000;           % Cooling water mass flow                      [kg / s]
C_cw   = c_cw * m_cw;       % Cooling water thermal capacity               [W / C]
h_cond = 500;               % Conduction coefficient                       [W / C]
h_conv = 0.1;               % Convection coefficient                       [W / C A] 
R_th   = 0.0817;            % Thermal resistance from [2]                  [K / W]
T_cwin = 20 + 0;            % Inlet cooling water temperature              [C]
V_tn   = 1.48;              % Thermo-neutral voltage                       [V]
tau_th = C_th * R_th;       % Thermal time constant                        [s]

%%                            HEAT RECOVERY

% COMBUSTION
% Post-combustion reaction rates
% H2 + 1/2 O2 ---> H2O           
eta_comb = 1;
xi4      = n_H2out * eta_comb;      % H2 comb. reaction rate               [mol / s]
Tcomb    = 590 + 273.15;            % Combustion temperature               [K]
pcomb    = p_an;                    % Combustion pressure                  [Pa]

% THERMAL LOAD
DeltaT_hx = 5;               % DeltaT of the hot water heat exchanger      [K]
T_HWin    = 45 + 273.15;     % Nominal inlet temp. of hot water from TU    [K]

% HEAT EXCHANGER
eta_hx = 0.92;               % Heat exchanger efficiency
tau_hx = 15.0;               % Heat exchanger time constant                [s]

% SOFC EXHAUST
Texh = T_HWin + DeltaT_hx;   % Exhaust temperature - depending on the load [K]
Tavg = (Tcomb+ Texh) / 2;

% Water quality
x_H2Ovap   = WaterSaturationPressure (Texh) / (pcomb/101325);              
h_H2Opcvap = WaterEnthalpyVaporization (pcomb/101325, Texh);             % [J/g]

% Enthalpies after combustion [J / kg]
C_H2comb  = refpropm('C', 'T', Tavg, 'P', pcomb/1000, 'HYDROGEN.FLD');
C_H2Ocomb = refpropm('C', 'T', Tavg, 'P', pcomb/1000, 'WATER.FLD');
h_aircomb = refpropm('H', 'T', Tcomb, 'P', pcomb/1000, 'AIR.PPF');
h_Arcomb  = refpropm('H', 'T', Tcomb, 'P', pcomb/1000, 'ARGON.FLD');
h_CO2comb = refpropm('H', 'T', Tcomb, 'P', pcomb/1000, 'CO2.FLD');
h_N2comb  = refpropm('H', 'T', Tcomb, 'P', pcomb/1000, 'NITROGEN.FLD');
h_O2comb  = refpropm('H', 'T', Tcomb, 'P', pcomb/1000, 'OXYGEN.FLD');

% Enthalpies at exhaust [J/kg]
h_airexh = refpropm('H', 'T', Texh, 'P', pcomb/1000, 'AIR.PPF');
h_Arexh  = refpropm('H', 'T', Texh, 'P', pcomb/1000, 'ARGON.FLD');
h_CO2exh = refpropm('H', 'T', Texh, 'P', pcomb/1000, 'CO2.FLD');
h_N2exh  = refpropm('H', 'T', Texh, 'P', pcomb/1000, 'NITROGEN.FLD');
h_O2exh  = refpropm('H', 'T', Texh, 'P', pcomb/1000, 'OXYGEN.FLD');

% Recovery enthalpies [J / mol] - post-combustion
h_airpc  = (h_aircomb - h_airexh) / 1000 * M_air; 
h_Arpc  = (h_Arcomb - h_Arexh) / 1000 * M_Ar; 
h_CO2pc = (h_CO2comb - h_CO2exh) / 1000 * M_CO2; 
h_H2pc  = C_H2comb * (Tcomb - Texh) / 1000 * M_H2; 
h_H2Opc = C_H2Ocomb * (Tcomb - Texh) / 1000 * M_H2O; 
h_N2pc  = (h_N2comb - h_N2exh) / 1000 * M_N2; 
h_O2pc  = (h_O2comb - h_O2exh) / 1000 * M_O2; 

%%                            LIBRARIES

%%                         IMPLEMENTED DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Relative humidity - function of inlet temperature
% c1    = [323.15, 353.15, 373.15]';
% c2    = [1, 0.50, 0.40]';
% polRH = fit(c1, c2, 'poly2');
% RH    = polRH(T);
% RH(T < 343.15) = 1;

% h_H2Orec_in = refpropm('H', 'T', T_ref, 'P', p_an * 100, 'WATER.FLD');
% h_H2O_T = refpropm('H', 'T', T, 'P', p_an * 100, 'WATER.FLD');

% lambda    = 0.90;   % Oxygen-to-hydrogen ratio; It is equal to 2 in [6]

% % Specific enthalpy for air pre-heating [J / mol]
% h_airph = C_air * (T - T_amb) * M_air / eta_ph;

% SPECIFIC VOLUME
% v_CH4   = 0.64828e03;     % Volumic mass of H2 - gas                       [g / m^3]
% v_H2    = 0.08988e03;     % Volumic mass of H2 - gas                       [g / m^3]
% v_H240  = 3.1777e03;      % Volumic mass of H2 at (p,T) = (40 bar, 25 C)   [g / m^3]
% v_H250  = 3.4443e03;      % Volumic mass of H2 at (p,T) = (50 bar, 70 C)   [g / m^3]
% v_H2O   = 999.97e03;      % Volumic mass of H2O - liquid                   [g / m^3]
% v_O2    = 1.308e03;       % Volumic mass of O2 - gas                       [g / m^3]
% v_O240  = 52.876e03;      % Volumic mass of H2 at (p,T) = (40 bar, 25 C)   [g / m^3]
% v_O250  = 56.761e03;      % Volumic mass of H2 at (p,T) = (50 bar, 70 C)   [g / m^3]

% disp('Hakuna Matata')