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

%%
function [V_, MF, VF, eta_e, P] = ALKEC_cellstack (A, I, T, T_amb, N_cell, ...
    p_an, p_cat, y_H2, y_H2O, y_O2, F, R, r1, r2, s, t1, t2, t3, f1, ...
    f2, C_th, h_cond, h_conv, R_th, U, theta)

%%                       ELECTROCHEMICAL MODEL

% PARTIAL PRESSURES
p_O2  = y_O2 *  p_an;    % O2 partial pressure                             [Pa]
p_H2  = y_H2 *  p_cat;   % H2 partial pressure                             [Pa]
p_H2O = y_H2O * p_an;    % H2O partial pressure                            [Pa]

% NERNST POTENTIAL [V]
E0 = 1.5184 - 1.5421e-03 * T + 9.523e-05 * T * log(T) + 9.84e-08 * T^2;   
                       % Ideal reversible potential at standard conditions [V]
                           
E = E0 + R * T * log(p_H2 * sqrt(p_O2) / p_H2O) / (2 * F); 

% REVERSIBLE VOLTAGE [V]
V_rev = E * theta;

% STACK VOLTAGE [V]
T_ = T - 273.15;                                       % Cell temperature  [C]
V_ = V_rev + (r1 + r2 * T_) * I / A + ...
    s * log10((t1 + t2 / T_ + t3 / T_^2) * I / A + 1);  
    

%%                         H2 PRODUCTION MODEL 

% FARADAY EFFICIENCY
etaF = (I / A)^2 / (f1 + (I / A)^2) * f2;

% FARADAY LAW - molar flows
n_H2  = etaF * N_cell * I / (2 * F);            % H2 molar flow            [mol / s]
n_H2O = 1 / U * N_cell * I * etaF / (2 * F);    % H2O molar flow           [mol / s]
n_O2  = N_cell * I * etaF / (4 * F);            % O2 molar flow            [mol / s]

% MASS FLOWS
M_H2   = 2.01588;                       % Molar mass of H2                 [g / mol] 
M_H2O  = 18.01528;                      % Molar mass of H2O                [g / mol]
M_O2   = 31.9988;                       % Molar mass of O2                 [g / mol]
MF_H2  = n_H2 * M_H2 / 1000;            % H2 mass flow                     [kg / s]
MF_H2O = n_H2O * M_H2O / 1000;          % H2O mass flow                    [kg / s]
MF_O2  = n_O2 * M_O2 / 1000;            % O2 mass flow                     [kg / s]
MF     = [MF_H2, MF_O2, MF_H2O];        % Mass flows                       [kg / s]

% VOLUMETRIC FLOWS
v_H2   = 0.08988e03;                    % Volumic mass of H2 - gas         [g / m^3]
v_H2O  = 999.97e03;                     % Volumic mass of H2O - liquid     [g / m^3]
v_O2   = 1.308e03;                      % Volumic mass of O2 - gas         [g / m^3]
VF_H2  = n_H2 * 3600 * M_H2 / v_H2;     % H2 volumetric flows              [m3 / h]
VF_H2O = n_H2O * 3600 *M_H2O / v_H2O;   % H2O volumetric flows             [m3 / h]
VF_O2  = n_O2 * 3600 * M_O2 / v_O2;     % O2 volumetric flows              [m3 / h]
VF     = [VF_H2, VF_O2, VF_H2O];        % Volumetric flows                 [m3 / h]

%%                     EFFICIENCY AND ELECTRICAL POWER

% THERMO-NEUTRAL VOLTAGE [V]
V_tn = 1.48;

% ELECTRICAL EFFICIENCY
eta_e = V_tn / V_;

% ABSORBED POWER 
eta_i = 0.98;                           % Inverter efficiency    
P     = N_cell * V_ * I;                % DC power                         [W]
P_AC  = P / eta_i;                      % AC power                         [W]

%%                            THERMAL MODEL

% THERMAL PARAMETERS
c_cw   = 4.18e03;                       % Cooling water specific heat      [J / kg C]
m_cw   = 0.0850;                        % Cooling water mass flow          [kg / s]
C_cw   = c_cw * m_cw;                   % Cooling water thermal capacity   [W / C]
UA     = h_cond + h_conv * I;           % Heat transfer coefficient        [W / C]
T_cwin = 20 + 0;                        % Inlet cooling water temperature  [C]
tau_th = C_th * R_th;                   % Thermal time constant            [s]

% THERMAL COEFFICIENTS
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