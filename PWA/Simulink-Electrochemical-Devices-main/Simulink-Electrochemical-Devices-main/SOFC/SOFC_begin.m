%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                                  SOFC                                   %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic modeling of a solid oxide fuel cell (SOFC). 
% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the solid oxide fuel cell are captured.

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

%%

clear variables
close all
clc;
tic;

%%                          INPUT PARAMETERS

% FLAGS
flag_fuel  = 1;              % 0: H2 as a fuel; 1: CH4 as a fuel
flag_dynam = 1;              % Dynamic calculation
flag_opcur = 1;              % Off-design calculation
flag_plot  = 1;              % Plot results
flag_appr  = 1;              % 0: no approximation, 1: PWA approximations of conversion curves

% MODELING PARAMETERS - Input stream
I_des  = 60;                 % Cell current at design                      [A]
F_f    = 1.0;                % Fuel factor: 1 for CH4, 4 for H2
n_fuel = 0.351e-02 * F_f;    % Fuel inlet flow at design conditions        [mol/s]
p_atm  = 1.01325e05;         % Atmosferic pressure                         [Pa]
p_an   = p_atm * 1.15;       % Anode pressure                              [Pa]
p_cat  = p_atm * 1.15;       % Cathode pressure                            [Pa] 
T      = 750 + 273.15;       % Cell temperature                            [K]
y_CH4  = 0.2832;             % CH4 molar fraction - anode inlet   
y_CO   = 0.0207;             % CO molar fraction - anode inlet
y_CO2  = 0;                  % CO2 molar fraction - anode inlet
y_C2H6 = 0;                  % C2H6 molar fraction - anode inlet
y_C3H8 = 0;                  % C3H8 molar fraction - anode inlet
y_H2   = 0.0482;             % H2 molar fraction - anode inlet
y_H2O  = 0.6348;             % H2O molar fraction - anode inlet
y_N2   = 0.0131;             % N2 molar fraction - cathode inlet
x_Ar   = 0.0092;             % Ar molar fraction - cathode inlet
x_CO2  = 0.0004;             % CO2 molar fraction - cathode inlet
x_H2O  = 0.0103;             % H2O molar fraction - cathode inlet
x_N2   = 0.7728;             % N2 molar fraction - cathode inlet
x_O2   = 0.2073;             % O2 molar fraction - cathode inlet
lambda = 9.15;               % Air-fuel factor 

% MODELING PARAMETERS - Cell geometry
N_cell = 38;                 % Number of cell in series
A_cell = 180;                % Cell active area                            [cm2]

% CONSTANTS
F     = 96485;               % Faraday's constant                          [C / mol]
R     = 8.314;               % Universal gas constant                      [J / mol-K]
T_a   = 25;                  % Ambient temperature                         [C]

%%                         SIMULATE TECHNOLOGY

n_fueldes = n_fuel;

if flag_fuel == 0
 
 n_H2in = n_fuel;
 run SOFC_H2data                      
 run SOFC_H2sim
 
else
    
 n_CH4in = n_fuel;
 run SOFC_CH4data                      
 run SOFC_CH4sim

end

%%                             PLOT RESULTS

if flag_plot == 1
 
 % ELECTRICAL AND THERMAL EFFICIENCIES
 figure;
 subplot(1,3,1)
 plot(P / 1000, etae_aux, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Electrical power output, {\it{P_{o}}} [kW]')
 ylabel('Electrical efficiency without auxiliaries, {\it{\eta_{e}}} [-]')
 subplot(1,3,2)
 plot(P / 1000, etat, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Thermal power input, {\it{Q_{i}}} [kW]')
 ylabel('Electrical power output with auxiliaries, {\it{P_{o}}} [kW]')
 subplot(1,3,3)
 plot(P./etae / 1000, P/1000, 'Linewidth', 1.5);
 set(gca, 'Fontsize', 12);
 xlabel('Thermal power input, {\it{Q_{i}}} [kW]')
 ylabel('Electrical power output without auxiliaries, {\it{P_{o}}} [kW]')
 
%  [p1,q1] = polyfit(P/1000./etae_aux, P/1000, 2);
%  r = polyval(p1, P/1000./etae_aux);
%  [p2,q2] = polyfit(P/1000./etae_aux, Q/1000, 2);
%  s = polyval(p2, P/1000./etae_aux);
%  p3 = [P(1)/1000./etae_aux(1), P(15)/1000./etae_aux(15), P(35)/1000./etae_aux(35), P(50)/1000./etae_aux(50)];
%  q3 = [P(1)/1000, P(15)/1000, P(35)/1000, P(50)/1000];
%  q4 = [Q(1)/1000, Q(15)/1000, Q(35)/1000, Q(50)/1000];

%  % EFFICIENCY
%  figure;
%  plot(P/1000./etae_aux, P/1000, P/1000./etae_aux, Q/1000,'Linewidth', 1.5);
%  hold on;
%  plot(P/1000./etae_aux, r, P/1000./etae_aux, s,'Linewidth', 1.5);
%  plot(p3, q3, p3, q4,'Linewidth', 1.5);
%  set(gca, 'Fontsize', 12);
%  xlabel('Inlet power, {\it{P_{in}}} [kW]')
%  ylabel('Outlet power, {\it{P_{out}}} [kW]')
 
%  % DYNAMIC
%  Pout_dyn.data(Pout_dyn.data < 0.0001) = 0.0001;
%  Qout_dyn.data(Qout_dyn.data < 0.0001) = 0.0001;
%  figure;
%  plot(t, Pout_dyn.data/1000, t-20, Qout_dyn.data/1000, 'Linewidth', 1.5);
%  set(gca, 'Fontsize', 12);
%  xlabel('Time, {\it{t}} [s]')
%  ylabel('Power output, {\it{P}} [kW]')
 
 % POLARIZATION CURVE
%  figure;
%  plot(J/A_cell, V/N_cell, 'Linewidth', 1.5);
%  xlabel('Current, {\it{I}} [A]');
%  ylabel('Voltage, {\it{V}} [V]');
%  set(gca, 'Fontsize', 12);
 
end

% REDUCED ORDER MODELS
if flag_appr == 1
  
  % Time constants [hr]
  Pdyn  = find(Pout_dyn.data >= 0.63 * Pout_dyn.data(end));
  tau_P = min(Pdyn) * 2 / 3600;
  Qdyn  = find(Qout_dyn.data >= 0.63 * Qout_dyn.data(end));
  tau_Q = min(Qdyn) * 2 / 3600;
  
  % PWA coefficients 1 kW
  N_bp     = 2;
  P_SOFC   = P / 1000 / 1.6;
  P_SOFCin = P_SOFC ./etae_aux;
  y1       = P_SOFCin(end);
  idx      = round(linspace(1, 50, N_bp));
  x_bp_val = P_SOFCin(idx)';
  y_bp_val = P_SOFC(idx)';
  figure;
  plot(P_SOFCin,P_SOFC,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  legend('Thermodynamic model', 'PWA approximation')
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,1) = m_etae;
  qq(:,1) = q_etae;
  
  p_SOFCin = P_SOFCin;
  p_SOFC   = P_SOFCin * mm + qq(1);
  
  MSE = sum((p_SOFC - P_SOFC).^2) / length(p_SOFC);
  
  % Power to heat ratio
  R      = P_SOFC ./ (P_SOFC ./ etae_aux .* etat);
  Q_SOFC = (P_SOFC ./ etae_aux .* etat);
  r     = polyfit([P_SOFC(1) P_SOFC(end)], [Q_SOFC(1), Q_SOFC(end)], 1);
  
  % PWA coefficients 100 kW
  P_SOFC   = P * 100 / 1000 / 1.6;
  P_SOFCin = P_SOFC ./ etae_aux;
  y2       = P_SOFCin(end);
  beta     = polyfit(P_SOFCin, P_SOFC, 2);
  x_bp_val = P_SOFCin(idx)';
  y_bp_val = P_SOFC(idx)';
  figure;
  plot(P_SOFCin,P_SOFC,'o', x_bp_val, y_bp_val, 'Linewidth', 1.5);
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
  Tech.SOFC.etal  = etae_aux(end);
  Tech.SOFC.alpha = polm;
  Tech.SOFC.beta  = polq(:,1);
  Tech.SOFC.taue  = tau_P;
  Tech.SOFC.taut  = tau_Q;
  Tech.SOFC.R     = 0.85;
  Tech.SOFC.aux   = P(1)/1000/1.6/etae_aux(1);
  save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
end

% CPU TIME
tCPU = toc;

%%                                     HAKUNA MATATA

% mm_(i,:)  = polyval(polm(i,:), [1 100]);
% qq_(i,:)  = polyval(polq(i,:), [1 100]);
% plot([1 100], mm_(i,:),'-', 'color', 'k');
% plot([1 100], mm(i,:),'--', 'color', 'k');
% plot([1 100], qq_(i,:),'-');
% plot([1 100], qq(i,:),'--');

% Tech.SOFC.pol_me = polm;
% Tech.SOFC.pol_qe = polq;
% Tech.SOFC.Lx     = x_bp_val(1:end-1);
% Tech.SOFC.Ux     = y_bp_val(2:end);
% Tech.SOFC.taue   = tau_P;
% Tech.SOFC.taut   = tau_Q;
% Tech.SOFC.R      = R;
% Tech.SOFC.etamax = max(etae_aux);
% Tech.SOFC.L      = 0.135;
