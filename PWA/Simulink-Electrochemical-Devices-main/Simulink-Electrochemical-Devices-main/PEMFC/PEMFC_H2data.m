%-------------------------------------------------------------------------%
%                    PEMFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                             Process Engineering Institute, September 2015


%-------------------------------------------------------------------------%
%                         PEMFC-CHP Parameters                            %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic modeling of a proton exchange membrane fuel cell (PEMFC). 
% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM fuel cell are captured.

% The modeled based on a generic example of PEM fuel cell using methane 
% (1.4 kW). The polarization curve is based on data from PSI cell and from 
% literature.
% The operative power can range in 0-120% of the rated power.
% The fuel can use either natural gas or hydrogen as a fuel. The operation
% pressure would be higher when using hydrogen, assuming to have it in a 
% storage at a pressure higher than the ambient pressure.
% A reference temperature of 80 °C is chosen as the commercial product is
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

% PRESSURES OF STORAGE
F_pst = 0;                % Storage pressure factor

%%                       THERMODYNAMIC DATA

% MOLAR WEIGHT
M_air = 28.97;            % air molar weight                               [kg / kmol]
M_H2  = 2;                % H2 molar weight                                [kg / kmol]
M_H2O = 18;               % H2O molar weight                               [kg / kmol]
M_CH4 = 16;               % CH4 molar weight                               [kg / kmol]
M_CO2 = 44;               % CO2 molar weight                               [kg / kmol]
M_N2  = 28;               % O2 molar weight                                [kg / kmol]
M_O2  = 32;               % N2 molar weight                                [kg / kmol]

% SPECIFIC HEAT
C_CH4    = 2.293;         % CH4 specific heat at (p, T) = (1 bar, 325 K)   [kJ / kg K]
C_CH4ref = 4.348;         % CH4 specific heat at (p, T) = (1 bar, 970 K)   [kJ / kg K]
C_H2     = 14.32;         % H2 specific heat at (p, T) = (1 bar, 350 K)    [kJ / kg K]
C_H2O    = 4.30;          % H2O specific heat at (p, T) = (1 bar, 350 K)   [kJ / kg K]
C_H2Opc  = 2.30;          % H2O specific heat at (p, T) = (1 bar, 700 K)   [kJ / kg K]
C_air    = 1.01;          % Air specific heat at (p, T) = (1 bar, 373 K)   [kJ / kg K]
C_airpc  = 1.10;          % Air specific heat at (p, T) = (1 bar, 700 K)   [kJ / kg K]

% LOWER HEATING VALUE
LHV_CH4   = 50e03;        % CH4 lower heating value                        [J/g]
LHV_H2    = 119.96e03;    % H2 lower heating value                         [J/g]

% SPECIFIC HEAT RATIO
gamma_air = 1.33;         % Air specific heat ratio
gamma_CH4 = 1.32;         % CH4 specific heat ratio
gamma_H2  = 1.409;        % H2 specific heat ratio

%%                           AUXILIARIES

% EFFICIENCIES
eta_ACDC = 0.92;          % Electrical efficiency of AC/DC converter
eta_ph   = 0.92;          % Pre-heater efficiency

% COMPRESSORS
beta_c = p_atm / p_an;    % Pressure ratio
eta_c  = 0.75;            % Compressor efficiency

% TIME CONSTANTS
tau_c  = 0.10;            % Compressor time constant                       [s]
tau_ph = 3.00;            % Air pre-heater time constant                   [s]

% MISCELLANEOUS AUXILIARIES
P_gamma = 155.0;          % Heat-up/Miscellaneous losses                   [W]

% CH4 specific compression work [J / mol]  
h_H2c  = gamma_H2 * R * T / (gamma_H2 - 1) * ((1 / beta_c)^...
    ((gamma_H2 - 1) / gamma_H2) - 1) / eta_c * F_pst;   

% Air specific compression work - air [J / mol]  
h_airc = gamma_air * R * T / (gamma_air - 1) * ((1 / beta_c)^...
    ((gamma_air - 1) / gamma_air) - 1) / eta_c; 

% Air pre-heating specific enthalpy [J / mol]
h_airph = C_air * (T - T_amb) * M_air / eta_ph;

% H2 pre-heating specific enthalpy [J / mol]
h_H2ph = C_H2 * (T - T_amb) * M_H2 / eta_ph;

% MOLAR FRACTIONS
y_H2O     = (1 - y_H2);      
y_CO2_ref = 1 - y_H2_ref - y_H2O_ref;
y_CO2     = 1 - y_H2 - y_H2O;

%%                           AIR HUMIDIFIER

% WATER SATURATION PRESSURE [atm]
pH2Osat = WaterSaturationPressure (T);

% AIR UTILIZATION - Data from Aspen simulations
split_air = 1.00;   % split factor for air: FC - Combustor

% RELATIVE HUMIDITY
RH = RelativeHumidity (T);

% WATER CONTENT
p_H2Ovap = RH * pH2Osat;      % Water partial pressure in feeding air     [atm]
x_H2O    = p_H2Ovap / p_cat;   % Water mole fraction in the feeding air
F_H2O    = 0.5;                % Water consumption factor

% WATER VAPORIZATION
p_WHSG    = p_an;              % WHSG pressure (flag when negative)        [atm]
T_H2Oinat = T_amb;             % Inlet water temperature air treatment     [K]
eta_WHSG  = 0.97;              % Thermal efficiency of WHSG
tau_WHSG  = 50.00;             % WHSG time constant                        [s]

% Water enthalpy of vaporization [J / g]
h_H2Ovap = WaterEnthalpyVaporization (p_WHSG, T) / (eta_WHSG * M_H2O);                     

% WH specific enthalpy [J / g]
h_H2OWH = C_H2O * (T - T_H2Oinat) / eta_WHSG;

% Overall specific enthalpy for WH + SG [J / g]
h_WHSG = h_H2OWH + h_H2Ovap;

%%                                  PEMFC STACK

% UTILIZATION FACTOR
c1   = [n_fueldes(1) * 0.1, n_fueldes(1), n_fueldes(1) * 1.25]';
c2   = [0.99, 0.90, 0.87]';
polU = polyfit(c1, c2, 2);
U    = polyval(polU, n_H2in(1));

% ANODE AND CATODE CHANNELS - Data from [3]
% Current
I = n_H2in(1) * 2 * F * U / N_cell;        % Operating current           [A]

% Molar flows after reforming - Data from Aspen simulations
n_H2Oin  = y_H2O * n_H2in(1) / y_H2;         % H2O inlet molar flow        [mol/s]
n_O2in   = lambda * n_H2in(1);               % O2 inlet molar flow         [mol/s]
n_airin  = n_O2in / (x_O2 * (1 - x_H2O));    % air inlet molar flow        [mol/s]
n_H2Ocin = n_airin * x_H2O;                  % water inlet molar flow      [mol/s]
n_H2r    = N_cell * I / (2 * F);             % H2 reacting molar flow      [mol/s]
n_H2Or   = N_cell * I / (2 * F);             % H2O reacting molar flow     [mol/s]
n_O2r    = N_cell * I / (4 * F);             % H2O reacting molar flow     [mol/s]
n_H2out  = n_H2in(1) - n_H2r;                % H2O outlet molar flow       [mol/s]
n_H2Oout = n_H2Oin + n_H2Or;                 % H2O outlet molar flow       [mol/s]
n_O2out  = n_O2in - n_O2r;                   % H2O outlet molar flow       [mol/s]

% Total mole flow at outlet
n_anout  = (n_H2Oout + n_H2out) / 0.85;      % Anode outlet molar flow     [mol/s]
n_catout = n_O2out + n_H2Ocin * F_H2O;       % Cathode outlet molar flow   [mol/s]

% Mole fractions at outlet
y_H2Oout = n_H2Oout / n_anout;               % H2O outlet molar fraction
y_H2out  = n_H2out / n_anout;                % H2 outlet molar fraction
x_O2out  = n_O2out / n_catout;               % O2 outlet molar fraction

% Partial pressure at inlet
p_H2in  = y_H2 *  p_an;                     % H2 inlet partial pressure    [Pa]
p_H2Oin = y_H2O *  p_an;                    % H2O inlet partial pressure    [Pa]
p_O2in  = x_O2 *  p_cat;                    % O2 inlet partial pressure    [Pa]

% Partial pressure at outlet
p_H2out  = y_H2out *  p_an;                 % H2 outlet partial pressure   [Pa]
p_H2Oout = y_H2Oout *  p_an;                % H2O outlet partial pressure  [Pa]
p_O2out  = x_O2out *  p_cat;                % O2 outlet partial pressure   [Pa]

% Valve constants
K_H2_  = n_H2out / p_H2out;                 % H2 valve constant            [mol/atm s]   
K_H2O_ = n_H2Oout / p_H2Oout;               % H2O valve constant           [mol/atm s]
K_O2_  = n_O2out / p_O2out;                 % O2 valve constant            [mol/atm s] 
K_H2   = 0.04220;                           % H2 valve constant            [mol/atm s]   
K_H2O  = 0.00772;                           % H2O valve constant           [mol/atm s]
K_O2   = 0.02110;                           % O2 valve constant            [mol/atm s] 

% Time constants
tau_H2  = 3.37;                             % H2 time constant             [s]
tau_H2O = 18.41;                            % H2O time constant            [s]
tau_O2  = 6.74;                             % O2 time constant             [s]

% Calculation parameters
K_r = N_cell / (4 * F);                     % Reaction parameter           [mol / C]

% NERNST VOLTAGE
E0  = 1.229 - 0.85e-03 * (T - 298.15);      % Ideal reversible potential   [V]   

% INITIAL PRESSURES
p_H20 = y_H2_ref * p_an;
p_O20 = x_O2 * p_cat;

% CONCENTRATION POLARIZATION LOSSES - Data from [2]
B     = 0.016;    % Fitting coefficient for concentration losses           [V]      
i_max = 1.2;      % Maximum current density                                [A / cm2] 

% ACTIVATION LOSSES - Data from [4]
% Concentration at membrane interface - Data from [1]
p_O2 = (p_O2in + p_O2out) / 2;
c_O2 = p_O2 / (5.08e06 * exp(-498 / T));   % O2 concentration at TPB       [mol / cm3]

% Fitting parameters - Data from [1]
xi1 = -0.949;
xi2 = 0.003420;
xi3 = 7.6e-05;
xi4 = -1.93e-04;

% OHMIC LOSSES - Data from [4]
R_el = 0.00016;             % Equivalent internal resistance               [ohm] 
z_m  = 178e-04;             % Membrane thickness                           [cm]
psi  = 23;   

% THERMAL MODEL
C_th   = 350e03;            % Thermal capacitance from [2]                 [J / K]
c_cw   = 4.18e03;           % Cooling water specific heat                  [J / kg C]
m_cw   = 0.25000;           % Cooling water mass flow                      [kg / s]
C_cw   = c_cw * m_cw;       % Cooling water thermal capacity               [W / C]
h_cond = 500;               % Conduction coefficient                       [W / C]
h_conv = 0.1;               % Convection coefficient                       [W / C A] 
R_th   = 0.0817;            % Thermal resistance from [2]                  [K / W]
T_cwin = 20 + 0;            % Inlet cooling water temperature              [C]
V_tn   = 1.48;              % Thermo-neutral voltage                       [V]
tau_th = C_th * R_th;       % Thermal time constant                        [s]


%%                                   HEAT RECOVERY

% COMBUSTION
Tcomb = 600.15;             % Combustion temperature                       [K]
pcomb = p_atm * 1.1;        % Combustion pressure                          [atm]

% THERMAL LOAD
DeltaT_hx = 5;              % DeltaT of the hot water heat exchanger       [K]
T_HWin    = 45 + 273.15;    % Nominal inlet temp. of hot water from TU     [K]

% HEAT EXCHANGER
eta_hx = 0.92;              % Heat exchanger efficiency
tau_hx = 15.0;              % Heat exchanger time constant                 [s]

% PEMFC EXHAUST
Texh = T_HWin + DeltaT_hx;  % Exhaust temperature - depending on the load  [K]

% Exhaust H2O flow [mol / s]
n_H2Oexh = n_H2Oin + n_H2Or + n_H2Ocin * F_H2O;
n_H2Oexh = n_H2Oexh + n_H2out;

% Exhaust N2 flow [mol / s]
n_N2exh  = n_O2in / x_O2 * (1 - x_O2);

% Exhaust flow [mol / s]
n_exh = (n_H2Oexh + n_N2exh) * (1 + 0.03);

% Water quality
x_H2Ocomb  = n_H2Oexh / n_exh; 
x_H2Ovap   = WaterSaturationPressure (Texh) / pcomb;
Delta_xH2O = x_H2Ocomb - x_H2Ovap;
n_H2Ovap   = Delta_xH2O * n_exh;
h_H2Opcvap = WaterEnthalpyVaporization (pcomb, T);

% Enthalpies after combustion [J / kg]
h_N2comb  = refpropm('H', 'T', Tcomb, 'P', pcomb * 100, 'NITROGEN.FLD');
h_O2comb  = refpropm('H', 'T', Tcomb, 'P', pcomb * 100, 'OXYGEN.FLD');
C_H2Ocomb = refpropm('C', 'T', Tcomb / 2, 'P', pcomb * 100, 'WATER.FLD');

% Exhaust enthalpies - water as a vapor [J / kg]
h_N2exh  = refpropm('H', 'T', Texh, 'P', pcomb * 100, 'NITROGEN.FLD');
h_O2exh  = refpropm('H', 'T', Texh, 'P', pcomb * 100, 'OXYGEN.FLD');

% Recovery enthalpies [J / mol] - post-combustion
h_N2pc  = (h_N2comb - h_N2exh) / 1000 * M_H2; 
h_H2Opc = C_H2Ocomb * (Tcomb - Texh) / 1000 * M_H2O; 
h_O2pc  = (h_O2comb - h_O2exh) / 1000 * M_O2; 

%%                            LIBRARIES

% WATER SATURATION PRESSURE
chi1 = 0.00129697;   % Fitting coefficient - Thermodynamic
chi2 = - 1.529053;   % Fitting coefficient - Thermodynamic
chi3 = 681.731481;   % Fitting coefficient - Thermodynamic
chi4 = - 136025.0;   % Fitting coefficient - Thermodynamic
chi5 = 10234070.7;   % Fitting coefficient - Thermodynamic

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

% disp('Hakuna Matata')