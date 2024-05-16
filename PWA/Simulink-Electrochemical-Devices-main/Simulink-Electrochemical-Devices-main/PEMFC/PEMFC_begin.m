%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                                 PEMFC                                   %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM fuel cell are captured.

% The modeled is based on a generic example of PEM fuel cell using methane 
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

%%

clear variables
close all
clc;
tic;

%%                                   INPUT PARAMETERS

% FLAGS
flag_fuel  = 0;              % 0: H2 as a fuel;     1: CH4 as a fuel
flag_oxi   = 1;              % 0: O2 as an oxidant; 1: air as an oxidant
flag_dynam = 1;              % Dynamic calculation
flag_opcur = 1;              % Off-design calculation
flag_plot  = 0;              % Plot results
flag_appr  = 2;              % 0: no approximation, 1: PWA approximations of conversion curves

% MODELING PARAMETERS - Input stream
I_des     = 30;              % Cell current at design                      [A]
n_fuel    = 2.51e-03;        % Fuel inlet flow at design conditions        [mol/s]
p_atm     = 1;               % Atmosferic pressure                         [atm]
p_an      = p_atm * 1.15;    % Anode pressure                              [atm]
p_cat     = p_atm * 1.50;    % Cathode pressure                            [atm]
pS_H2     = 40;              % H2 storage pressure                         [atm]
pS_O2     = 40;              % O2 storage pressure                         [atm]
T         = 70 + 273.15;     % Cell temperature                            [K]
T_a       = 25;              % Ambient temperature                         [C]
TS_H2     = 25 + 273.15;     % H2 storage temperature                      [K]
TS_O2     = 25 + 273.15;     % O2 storage temperature                      [K]
y_H2      = 0.995;           % H2 inlet molar fraction
y_H2_ref  = 0.66;            % H2 inlet molar fraction - reformer
y_H2O_ref = 0.15;            % H2O inlet molar fraction - reformer
x_O2      = 0.21;            % O2 inlet molar fraction
lambda    = 0.90;            % O2-to-H2 ratio; it is equal to 2 in [6]

% MODELING PARAMETERS - Cell geometry
N_cell = 72;                 % Number of cell in series
A_cell = 180;                % Cell active area                            [cm2]

% CONSTANTS
F     = 96485;               % Faraday's constant                          [C / mol]
R     = 8.314;               % Universal gas constant                      [J / mol-K]

%%                         SIMULATE TECHNOLOGY

n_fueldes = n_fuel;

if flag_fuel == 0
 
  if flag_oxi == 1
    n_H2in = n_fuel;
    run PEMFC_H2data                      
    run PEMFC_H2sim
  else
    % INPUT DATA FOR MODEL VALIDATION - Belenos/PSI B25 
    A_cell    = 238;            % cell active area       [cm2]
    I_des     = 0.8;              % design current density [A/cm2]
    N_cell    = 178;            % number of cells    
    p_atm     = 1;              % atmosferic pressure    [atm]
    p_an      = p_atm * 3.50;   % anode pressure         [atm]
    p_cat     = p_atm * 3.50;   % cathode pressure       [atm] 
    rH_H2     = 0.54;           % H2 relative humidity
    rH_O2     = 0.54;           % O2 relative humidity
    T         = 74 + 273.15;    % cell temperature       [K]
    z_m       = 30e-04;         % membrane dry thickness [cm]
    lambda_H2 = 1.3;
    lambda_O2 = 1.5;
    
    % SIMULATION DATA - System BoP
    I = I_des;
    run PEMFC_H2O2data
    
    % SIMULATION
    II = linspace(0.0025, 0.9, 50);
    for i = 1 : length(II)
      I = II(i);
      run PEMFC_H2O2data
      t = 0 : 0.1 : 120;
      nH2st(1:length(t)) = nH2r;
      nO2st(1:length(t)) = nO2r;
      sim('PEMFC_smlkH2O2', t, [], [t; [nH2st; nO2st]]')
      P(i)         = Psystem.data(end);
      p(i)         = Pstack.data(end);
      V(i)         = Vstack.data(end) / N_cell;
      Q(i)         = Qsystem.data(end);
      H(i)         = nH2r * M_H2 * LHV_H2;
      etastack(i)  = p(i) / H(i);
      etae_aux(i)   = P(i) / H(i) * 0.9;
      etae(i)      = etae_aux(i);
      P_(i)        = Psystem.data(end) / 25;
      H_(i)        = P_(i) / etae_aux(i);

    end
    
      plot(H_/1000,P_/1000)
      grid on;
      box on;
      figure;
      plot(P/1000,etae_aux)
      
  end
  
else
    
 n_CH4in = n_fuel;
 run PEMFC_CH4data                      
 run PEMFC_CH4sim

end
      
%%                             PLOT RESULTS

if flag_plot == 1
 
 % ELECTRICAL AND THERMAL EFFICIENCIES
 figure;
 subplot(1,3,1)
 plot(P / 1000, etae_aux, P / 1000, etat,'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Electrical power output, {\it{P_{o}}} [kW]')
 ylabel('Electrical efficiency, {\it{\eta_{e}}} [-]')
 subplot(1,3,2)
 plot(P./etae_aux / 1000, P/1000, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Thermal power input, {\it{Q_{i}}} [kW]')
 ylabel('Electrical power output with auxiliaries, {\it{P_{o}}} [kW]')
 subplot(1,3,3)
 plot(P./etae / 1000, P/1000, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Thermal power input, {\it{Q_{i}}} [kW]')
 ylabel('Electrical power output without auxiliaries, {\it{P_{o}}} [kW]')
 
%  figure;
%  plot(nfuel/n_fueldes, etae, nfuel/n_fueldes, etat, 'Linewidth', 1.5);
%  set(gca, 'Fontsize', 12);
%  xlabel('Normalized inlet fuel, {\it{n_f}} [-]')
%  ylabel('Conversion efficiency, {\it{\eta}} [-]')
%  legend('Electrical efficiency {\it{\eta_e}}','Thermal efficiency {\it{\eta_t}}')

%  % POLARIZATION CURVE
%  figure;
%  plot(II/A_cell/1e04, V/N_cell, 'Linewidth', 1.5);
%  set(gca, 'Fontsize', 12);
%  xlabel('Current density, {\it{j}} [A/m^2]')
%  ylabel('Cell voltage, {\it{V}} [V]')
 
end

% REDUCED ORDER MODELS
if flag_appr == 1
  
  % Time constants [hr]
%   Pdyn  = find(Pout_dyn.data >= 0.63 * Pout_dyn.data(end));
%   tau_P = min(Pdyn) * 2 / 3600;
%   Qdyn  = find(Qout_dyn.data >= 0.63 * Qout_dyn.data(end));
%   tau_Q = min(Qdyn) * 2 / 3600;
  
  % PWA coefficients 1 kW
  N_bp      = 2;
  P_PEMFC   = P / 1000 / 1.28;
  P_PEMFCin = P_PEMFC ./ etae_aux * 0.95;
  y1        = P_PEMFCin(end);
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMFCin(idx)';
  y_bp_val  = P_PEMFC(idx)';
  plot(P_PEMFCin,P_PEMFC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,1) = m_etae;
  qq(:,1) = q_etae;
  
  p_PEMFCin = P_PEMFCin;
  p_PEMFC   = P_PEMFCin * mm + qq(1);
  
  MSE = sum((p_PEMFC - P_PEMFC).^2) / length(p_PEMFC);
  
  % Power to heat ratio
  R       = P_PEMFC ./ (P_PEMFC ./ etae_aux .* etat);
  Q_PEMFC = (P_PEMFC ./ etae_aux .* etat);
  r       = polyfit([P_PEMFC(1) P_PEMFC(end)], [Q_PEMFC(1), Q_PEMFC(end)], 1);
  plot(P_PEMFC,Q_PEMFC);
  
  % PWA coefficients 100 kW
  P_PEMFC   = P * 100 / 1000 / 1.28;
  P_PEMFCin = P_PEMFC ./ etae_aux * 0.95;
  y2        = P_PEMFCin(end);
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMFCin(idx)';
  y_bp_val  = P_PEMFC(idx)';
  plot(P_PEMFCin,P_PEMFC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,2) = m_etae;
  qq(:,2) = q_etae;
  
  polm = zeros(size(mm,1),1);
  polq = zeros(size(mm,1),2);
    
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit([y1 y2], mm(i,:), 0);
    polq(i,:) = polyfit([y1 y2], qq(i,:), 1);
  end

  % SAVING RESULTS
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
  Tech.PEMFC.etal  = etae_aux(end) * 0.95;
  Tech.PEMFC.alpha = polm;
  Tech.PEMFC.beta  = polq(:,1);
  Tech.PEMFC.taue  = tau_P;
  Tech.PEMFC.taut  = tau_Q;
  Tech.PEMFC.R     = 0.95;
  Tech.PEMFC.aux   = P(1)/1000/etae_aux(1);
  save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
  
elseif flag_appr == 2
  
  % PWA coefficients 1 kW
  N_bp      = 2;
  P_PEMFC   = P_ / 1000;
  P_PEMFCin = H_ / 1000;
%   P_PEMFCin = H_ / 1000 / 1.15;
  y1        = P_PEMFCin(end);
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMFCin(idx)';
  y_bp_val  = P_PEMFC(idx)';
  plot(P_PEMFCin,P_PEMFC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,1) = m_etae;
  qq(:,1) = q_etae;
  
    % PWA coefficients 100 kW
  N_bp      = 5;
  P_PEMFC   = P_ / 1000 * 100;
  P_PEMFCin = H_ / 1000 * 100;
%   P_PEMFCin = H_ / 1000 / 1.15;
  y2        = P_PEMFCin(end);
  idx       = round(linspace(1, 50, N_bp));
  x_bp_val  = P_PEMFCin(idx)';
  y_bp_val  = P_PEMFC(idx)';
  plot(P_PEMFCin,P_PEMFC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,2) = m_etae;
  qq(:,2) = q_etae;
  
  polm = zeros(size(mm,1),1);
  polq = zeros(size(mm,1),2);
    
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit([y1 y2], mm(i,:), 0);
    polq(i,:) = polyfit([y1 y2], qq(i,:), 1);
  end
  
  % SAVING RESULTS H2
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
  Tech.H2PEMFC.etal  = P_(end) / H_(end);
  Tech.H2PEMFC.alpha = polm;
  Tech.H2PEMFC.beta  = polq(:,1);
  Tech.H2PEMFC.taue  = 2 / 3600;
  Tech.H2PEMFC.taut  = 25 / 3600;
  Tech.H2PEMFC.R     = 0.95;
  Tech.H2PEMFC.aux   = H_(1) / 1000;
  save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
  
%   % SAVING RESULTS O2
%   load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
%   Tech.O2PEMFC.etal  = P_(end) / H_(end);
%   Tech.O2PEMFC.alpha = polm;
%   Tech.O2PEMFC.beta  = polq(:,1);
%   Tech.O2PEMFC.taue  = 2 / 3600;
%   Tech.O2PEMFC.taut  = 25 / 3600;
%   Tech.O2PEMFC.R     = 0.95;
%   Tech.O2PEMFC.aux   = H_(1) / 1000;
%   save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
   
end

% CPU TIME
tCPU = toc;

%%                                    HAKUNA MATATA

%   Tech.PEMFC.pol_me = polm;
%   Tech.PEMFC.pol_qe = polq;
%   Tech.PEMFC.Lx     = x_bp_val(1:end-1);
%   Tech.PEMFC.Ux     = y_bp_val(2:end);
%   Tech.PEMFC.taue   = tau_P;
%   Tech.PEMFC.taut   = tau_Q;
%   Tech.PEMFC.R      = R;
%   Tech.PEMFC.etamax = max(etae_aux);
%   Tech.PEMFC.L      = 0.241;