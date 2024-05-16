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

%%                        THERMODYNAMIC MODELS

% STEADY-STATE POLARIZATION CURVE
if flag_dynam == 0      
 
 [V_, MF, VF, eta_e, P] = PEMEC_cellstack (I, T, p_an, p_cat, ...
    y_H2, y_O2, F, R, CH, cH2_me0, cO2_me0, DH, epsilon, ...
    epsilon_p, io_an, io_cat, Req_el, t_an, t_cat, t_mem, A, N_cell, ...
    Epro, T_ref, T_amb, C_th, R_th, h_cond, h_conv, U, p_H2, p_O2, p_H2O);

% DESIGN CALCULATION
elseif flag_opcur == 0
 
 % Time horizon
 t  = 0 : 1 : 500;
 z  = length(t);
 tf = t(z);
 
 % Inlet fuel time profile
 I(1:z) = I_des;
 
 % Simulation
 sim('PEMEC_smlk', t, [], [t; I]');
 
 % Design hydrogen/oxygen production
 VH2 = VF_H2.Data(end);
 VO2 = VF_O2.Data(end);
 P   = PP.Data(end);
 V   = VV.Data(end);
 
% OFF-DESIGN CALCULATION
elseif flag_opcur == 1
 
 % Time horizon
 t  = 0 : 1 : 500;
 z  = length(t);
 tf = t(z);
  
 % Pre-allocation for speeding
 II  = linspace(I_des * 0.05, I_des, 50);
 w   = length(II);
 P   = zeros(1, w);
 V   = zeros(1, w);
 VH2 = zeros(1, w);
 VO2 = zeros(1, w);

 for j = 1 : w
   
  % Inlet fuel time profile
  I(1:z) = II(j);
   
  % Simulation
  run PEMEC_data
  sim('PEMEC_smlk', t, [], [t; I]');
 
  % Off-design hydrogen/oxygen production
  QH2(j)      = Q_H2.Data(end);
  VH2(j)      = VF_H2.Data(end);
  VO2(j)      = VF_O2.Data(end);
  P(j)        = PP.Data(end);
  V(j)        = VV.Data(end);
  etaH(j)     = QH2(j) / ((P(j) - P_gamma) / 1000);
  etaH_aux(j) = QH2(j) / (P(j) / 1000);
  
 end
  
end