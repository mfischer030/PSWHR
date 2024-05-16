%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                      PEMEC - Polarization curve                         %
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

%%
function  [V_, MF, VF, eta_e, P] = PEMEC_cellstack (I, T, p_an, p_cat, ...
    y_H2, y_O2, F, R, CH, cH2_me0, cO2_me0, DH, epsilon, ...
    epsilon_p, io_an, io_cat, Req_el, t_an, t_cat, t_mem, A, N_cell, ...
    Epro, T_ref, T_amb, C_th, R_th, h_cond, h_conv, U, p_H2, p_O2, p_H2O)


%%                       ELECTROCHEMICAL MODEL

% TEMPERATURES
T_an  = T;                % Anode temperature                              [K]               
T_cat = T;                % Cathode temperature                            [K] 

% PARTIAL PRESSURES
p_O2  = p_O2 * 101325;    % O2 partial pressure                            [Pa]
p_H2  = p_H2 * 101325;    % H2 partial pressure                            [Pa]
p_H2O = p_H2O * 101325;   % H2O partial pressure                           [Pa]

% NERNST POTENTIAL [V]
E0 = 1.5184 - 1.5421e-03 * T + 9.523e-05 * T * log(T) + 9.84e-08 * T^2;   
                       % Ideal reversible potential at standard conditions [V]
                           
E = E0 + R * T * log(p_H2 * sqrt(p_O2) / p_H2O) / (2 * F); 

% DIFFUSION OVERPOTENTIAL [V]
a = 3.64e-04;   % Coefficient from [3]
b = 2.334;      % Coefficient from [3]

% H2-H2O diffusion coefficient [cm2 / s]
D_H2 = a * ((T_cat / sqrt(12.8 * 218.3)) ^ b) * ((33.3 * 647.3)^(1 / 3))...
    * ((12.8 * 218.3) ^ (5/12)) * ((1/2 + 1/18) ^ (1/2)) / (p_cat*1e-05);  

% O2-H2O diffusion coefficient [cm2 / s]
D_O2 = a * ((T_an / sqrt(49.7 * 218.3)) ^ b) * ((154.4 * 647.3) ^ (1/3))...
    * ((49.7 * 218.3) ^ (5/12)) * ((1/32 + 1/18) ^ (1/2)) / (p_an*1e-05);  

% H2-H2O effective diffusion coefficient  [cm2 / s]
Deff_cat = D_H2 * epsilon * (((epsilon - epsilon_p) / ...
    (1 - epsilon_p)) ^ 0.785); 

% O2-H2O effective diffusion coefficient  [cm2 / s]
Deff_an  = D_O2 * epsilon * (((epsilon - epsilon_p) / ...
    (1 - epsilon_p)) ^ 0.785);   
  
% Unit conversion - cm2 ---> m2
Deff_an  = Deff_an  * 1e04;   
Deff_cat = Deff_cat * 1e04; 

% Concentration at membrane interface [mol / m3]
CH2_me = p_cat * y_H2 / (R * T_cat) + t_cat * I / (2 * F * A * Deff_cat);
CO2_me = p_an * y_O2 / (R * T_an) + t_an * I / (4 * F * A * Deff_an);

% Overvoltage [V]
V_diff_an  = R * T_an / (4 * F) * log(CO2_me / cO2_me0);
V_diff_cat = R * T_cat / (4 * F) * log(CH2_me / cH2_me0);
V_diff     = V_diff_an + V_diff_cat;

% ACTIVATION OVERPOTENTIAL [V]
alpha_an  = 2;     % Anode charge tranfer coefficient
alpha_cat = 0.5;   % Cathode charge tranfer coefficient

% Calculation parameters
K1 = R * T_an / (alpha_an * F);
K2 = 1 / (2 * A * io_an);
K3 = R * T_cat / (alpha_cat * F);
K4 = 1 / (2 * A * io_cat);

% Overvoltage [V]
V_act_an  = K1 * asinh(K2 * I);
V_act_cat = K3 * asinh(K4 * I);
V_act     = V_act_an + V_act_cat;

% OHMIC OVERPOTENTIAL [V]
% Membrane conductivity [S / m]
sigma_mem  = 16 * exp(-Epro / R * (1 / T - 1 / T_ref));
sigma_mem1 = F^2 * CH * DH / (R * T);

% Overvoltage [V]
V_ohm = (Req_el + t_mem / (A * sigma_mem)) * I;

% TOTAL VOLTAGE
V_ = E + V_diff + V_act + V_ohm;   % Cell voltage                          [V]
V  = V_ * N_cell;                  % Stack voltage                         [V]

% ABSORBED POWER [W]
P = V * I;

% ELECTRICAL EFFICIENCY
V_tn  = 1.48;                      % Thermo-neutral efficiency             [V]
eta_e = V_tn / V_;                 % Electrical efficiency

%%                         H2 PRODUCTION MODEL 

etaF = 0.98;   % Faraday efficiency - value from Valverde et al. (2011)

% FARADAY LAW - molar flows
n_H2  = N_cell * I * etaF / (2 * F);     % H2 molar flow produced          [mol / s]
n_O2  = N_cell * I * etaF / (4 * F);     % O2 molar flow produced          [mol / s]
n_H2O = N_cell * I * etaF / (2 * F);     % H2O molar flow consumed         [mol / s]

% WATER TRANSPORT
% Constant parameters - Values from [3]
Dw      = 1.28e-10;                      % Water diffusion coefficient     [m2 / s]
Kw      = 1.58e-18;                      % Darcy's constant                [m2]
nd      = 7;                             % Electro-osmotic drag coeff.     [molH2O / molH+]
M_H2    = 2.01588;                       % Molar mass of H2                [g / mol] 
M_H2O   = 18.01528;                      % Molar mass of H2O               [g / mol]
M_O2    = 31.9988;                       % Molar mass of O2                [g / mol]
v_H2    = 0.08988e03;                    % Volumic mass of H2 - gas        [g / m^3]
v_H2O   = 999.97e03;                     % Volumic mass of H2O - liquid    [g / m^3]
v_O2    = 1.308e03;                      % Volumic mass of O2 - gas        [g / m^3]
rho_H2O = 1000;                          % Water density                   [kg / m3]
mu_H2O  = 1.1e-03;                       % Water dynamic viscosity         [Pa s]

% Water transport for electro-osmotic drag [mol / s]
n_H2Oeo = nd * I / F;  

% Water transport for pressure gradient [mol / s]
n_H2Opg = Kw * A * rho_H2O / (mu_H2O * M_H2O * t_mem) * (p_cat - p_an);

% Overall water through the membrane [mol / s]
n_H2Omem1 = n_H2Oeo - n_H2Opg + ...
    A * Dw / t_mem * (t_an * I / (2 * F * A * Deff_an));
n_H2Omem2 = 1 - Dw / t_mem * (t_cat / Deff_cat + t_an / Deff_an);
n_H2Omem  = n_H2Omem1 / n_H2Omem2;

% Water balances
n_H2Oan  = (n_H2O + n_H2Omem) * 1 / U;
n_H2Ocat = n_H2Omem; 

% MASS FLOWS
MF_H2     = n_H2 * M_H2 / 1000;                    % H2 mass flow          [kg / s]
MF_O2     = n_O2 * M_O2 / 1000;                    % O2 mass flow          [kg / s]
MF_H2Oan  = n_H2Oan * M_H2O / 1000;                % H2O mass flow - anode [kg / s]
MF_H2Ocat = n_H2Ocat * M_H2O / 1000;               % H2O mass flow - cat.  [kg / s]
MF        = [MF_H2, MF_O2, MF_H2Oan, MF_H2Ocat];   % Mass flows            [kg / s]

% VOLUMETRIC FLOWS
VF_H2     = n_H2 * 3600 * M_H2 / v_H2;             % H2 volumetric flow    [m3 / h] 
VF_O2     = n_O2 * 3600 * M_O2 / v_O2;             % O2 volumetric flow    [m3 / h]
VF_H2Oan  = n_H2Oan * 3600 *M_H2O / v_H2O;         % H2O vol. flow - an.   [m3 / h] 
VF_H2Ocat = n_H2Ocat * 3600 *M_H2O / v_H2O;        % H2O vol. flow - cat.  [m3 / h] 
VF        = [VF_H2, VF_O2, VF_H2Oan, VF_H2Ocat];   % Volume flows          [m3 / h]

%%                            THERMAL MODEL

% THERMAL PARAMETERS
c_cw   = 4.18e03;                       % Cooling water specific heat      [J / kg C]
m_cw   = 0.1075;                        % Cooling water mass flow          [kg / s]
C_cw   = c_cw * m_cw;                   % Cooling water thermal capacity   [W / C]
UA     = h_cond + h_conv * I;           % Heat transfer coefficient        [W / C]
T_cwin = 20 + 0;                        % Inlet cooling water temperature  [C]
tau_th = C_th * R_th;                   % Thermal time constant            [s]

% THERMAL COEFFICIENTS
%--------------------------------------------------------------------------
% T_      = T - 273.15;
% Q_gen   = N_cell * (V_ - V_tn) * I;
% Q_loss  = 1 / R_th * (T_ - T_amb);
% T_cwout = T_cwin + (T_ - T_cwin) * (1 - exp(-UA / C_cw));
% LMTD    = ((T_ - T_cwin) - (T_ - T_cwout)) / ...
%     (log((T_ - T_cwin) / (T_ - T_cwout)));
% Q_cool  = C_cw * (T_cwout - T_cwin);
% Q_tot   = Q_gen - Q_loss - Q_cool;
% dTdt    = Q_tot / C_th;
%--------------------------------------------------------------------------

a = 1 / tau_th + C_cw / C_th * (1 - exp(-UA / C_cw));
b = N_cell * V_ * I * (1 - eta_e) / C_th + T_amb / tau_th + ...
    C_cw * T_cwin / C_th * (1 - exp(-UA / C_cw));

TT = T - 273.15;

% TEMPERATURE VARIATION - slope [K / s]
dT_dt = -a * TT + b;


%--------------------------------------------------------------------------
% TAFEL'S APPROXIMATION
% if flag_Tafel == 0
%     
%  V_act_an  = K1 * asinh(K2 * I);
%  V_act_cat = K3 * asinh(K4 * I);
% 
% elseif flag_Tafel == 1
% 
%  % Tafel approximation
%  V_act_an  = K1 * log(K2 * I);  
%  V_act_cat = K3 * log(K4 * I);  
%  
% end
%--------------------------------------------------------------------------

end