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

%%

clear variables
close all
clc;
tic;

%%                           INPUT PARAMETERS

% FLAGS
flag_dynam = 1;         % Dynamic calculation
flag_opcur = 1;         % Off-design calculation
flag_plot  = 1;         % Plot results
flag_appr  = 1;         % 0: no approximation, 1: PWA approximations of conversion curves

% MODELING PARAMETERS - Input stream
I_des = 240;            % Operating current;                               [A]
p_atm = 101325;         % Atmosferic pressure                              [Pa]
p_an  = p_atm * 50;     % Anode pressure                                   [Pa]
p_cat = p_atm * 50;     % Cathode pressure                                 [Pa] 
T     = 80 + 273.15;    % Cell temperature                                 [K]
T_amb = 25;             % Ambient temperature                              [C]  
T_ref = 327.15;         % Reference temperature                            [K]
y_O2  = 0.005;          % O2 molar fraction                  
y_H2  = 0.005;          % H2 molar fraction                
y_H2O = 1 - y_O2;       % H2O molar fraction 
              
% MODELING PARAMETERS - Cell geometry
A_cell = 0.0170;        % MEA area                                         [m2]
N_cell = 200;           % Number of cell in series                   

% CONSTANT
F   = 96485;            % Faraday's constant                               [C/mol]
R   = 8.314;            % Universal gas constant                           [J/mol K]
T_a = 25;               % Ambient temperature                              [C]

%%                         SIMULATE TECHNOLOGY

I = I_des;
run PEMEC_data
run PEMEC_sim

%%                             PLOT RESULTS

if flag_plot == 1
 
%  % PRODUCED HYDROGEN AND OXYGEN
 hold on;
 plot(P / 1000, VH2, P/1000, VO2,  'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Electrical power input, {\it{P_{i}}} [kW]')
 ylabel('Gas production, {\it{V_j}} [m^3/hr]')
 legend('Hydrogen production', 'Oxygen production', 'Location','northwest')
 
 figure;
 subplot(1,2,1)
 plot(P / 1000, QH2,'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Electrical power input, {\it{P_{i}}} [kW]')
 ylabel('Hydrogen output, {\it{Q_{H_2}}} [-]')
 subplot(1,2,2)
 plot(QH2 ./ etaH, QH2, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Electrical power input, {\it{P_{i}}} [kW]')
 ylabel('Hydrogen output, {\it{Q_{H_2}}} [-]')
 
 figure;
 plot(P/1e05, VH2, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Normalized inlet power, {\it{P_e}} [-]')
 ylabel('Produced hydrogen, {\it{V_{H_2}}} [-]')

 % POLARIZATION CURVE
 figure;
 plot(II/A_cell/1e04, V/N_cell, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Current density, {\it{j}} [A/m^2]')
 ylabel('Cell voltage, {\it{V}} [V]')

end

% REDUCED ORDER MODELS
if flag_appr == 1
  
  % Time constants [hr]
  Hdyn  = find(n_H2.data >= 0.63 * n_H2.data(end));
  tau_H = min(Hdyn) * 2 / 3600;
  
  % PWA coefficients 100 kW
  N_bp      = 8;
  QH2_PEMEC = QH2;
%   QH2_PEMEC = QH2 * 0.5;
  P_PEMEC   = QH2 ./ etaH_aux;
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMEC(idx)';
  y_bp_val  = QH2_PEMEC(idx)';
  plot(P_PEMEC, QH2_PEMEC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,2) = m_etae;
  qq(:,2) = q_etae;
  
  % PWA coefficients 1 kW
  QH2_PEMEC = QH2 / 100;
%   QH2_PEMEC = QH2 / 100 * 0.5;
  P_PEMEC   = QH2 / 100 ./ etaH_aux;
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMEC(idx)';
  y_bp_val  = QH2_PEMEC(idx)';
  plot(P_PEMEC, QH2_PEMEC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,1) = m_etae;
  qq(:,1) = q_etae;
  
  polm = zeros(size(mm,1),1);
  polq = zeros(size(mm,1),2);
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit([1 100], mm(i,:), 0);
    polq(i,:) = polyfit([1 100], qq(i,:), 1);
  end

  %SAVING RESULTS
  %load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
  %load('C:\Users\gahu\Desktop\EMPA_gh\8_Coding');
  %Tech.H2PEMEC.etal  = etaH_aux(end);
  %Tech.H2PEMEC.alpha = polm;
  %Tech.H2PEMEC.beta  = polq(:,1);
  %Tech.H2PEMEC.taue  = tau_H;
  %Tech.H2PEMEC.taut  = [];
  %Tech.H2PEMEC.R     = [];
  %Tech.H2PEMEC.aux   = QH2(1) / 100 / etaH_aux(1);
  %save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
  %save('C:\Users\gahu\Desktop\EMPA_gh\8_Coding');
  
  % SAVING RESULTS
  %load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
  %Tech.O2PEMEC.etal  = etaH_aux(end);
  %Tech.O2PEMEC.alpha = polm;
  %Tech.O2PEMEC.beta  = polq(:,1);
  %Tech.O2PEMEC.taue  = tau_H;
  %Tech.O2PEMEC.taut  = [];
  %Tech.O2PEMEC.R     = [];
  %Tech.O2PEMEC.aux   = QH2(1) / 100 / etaH_aux(1);
  %save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
  
end

% CPU TIME
tCPU = toc;

%%                                    HAKUNA MATATA

% Tech.PEMEC.pol_me = polm;
% Tech.PEMEC.pol_qe = polq;
% Tech.PEMEC.Lx     = x_bp_val(1:end-1);
% Tech.PEMEC.Ux     = y_bp_val(2:end);
% Tech.PEMEC.taue   = tau_H;
% Tech.PEMEC.taut   = [];
% Tech.PEMEC.R      = [];
% Tech.PEMEC.etamax = [];
% Tech.PEMEC.L      = 0.0653;