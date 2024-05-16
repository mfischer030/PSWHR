%--------------------------------------------------------------------------------------------------%
%                                      H2-O2 PEM Fuel Cell                                         %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, July 2016

%--------------------------------------------------------------------------------------------------%
%                                 Cell stack fitting coefficients                                  %
%--------------------------------------------------------------------------------------------------%

clear variables
close all
clc

flag_plot_pol = 0;

%%                                 LOAD POLARIZATION CURVES

%load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\PEMFC\B25_T25.mat');
%load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\PEMFC\B25_T35.mat');
%load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\PEMFC\B25_T55.mat');
%load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\PEMFC\B25_T75.mat');
%load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\PEMFC\B25_T95.mat');

load('C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\ressources\Simulink-Electrochemical-Devices-main\Simulink-Electrochemical-Devices-main\PEMFC\B25_T25.mat');
load('C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\ressources\Simulink-Electrochemical-Devices-main\Simulink-Electrochemical-Devices-main\PEMFC\B25_T35.mat');
load('C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\ressources\Simulink-Electrochemical-Devices-main\Simulink-Electrochemical-Devices-main\PEMFC\B25_T55.mat');
load('C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\ressources\Simulink-Electrochemical-Devices-main\Simulink-Electrochemical-Devices-main\PEMFC\B25_T75.mat');
load('C:\Users\fism\Desktop\MA_thesis\02_modeling_and_optimization\ressources\Simulink-Electrochemical-Devices-main\Simulink-Electrochemical-Devices-main\PEMFC\B25_T95.mat');

%%                                   PLOT POLARIZATION CURVES

% COMPUTE ORIGINAL DATA
I.T25 = B25_T25(:,1);
I.T35 = B25_T35(:,1);
I.T55 = B25_T55(:,1);
I.T75 = B25_T75(:,1);
I.T95 = B25_T95(:,1);
V.T25 = B25_T25(:,2);
V.T35 = B25_T35(:,2);
V.T55 = B25_T55(:,2);
V.T75 = B25_T75(:,2);
V.T95 = B25_T95(:,2);
T     = [25 35 55 75 95] + 273.15;

if flag_plot_pol  == 1
  
  for i = 1 : length(T)

    a1 = fieldnames(I);
    a2 = char(a1(i));
    J  = extractfield(I, a2)';
    b1 = fieldnames(V);
    b2 = char(b1(i));
    U  = extractfield(V, b2)'; 
    
    hold on;
    plot(J, U, 'Linewidth', 2)
    set(gca, 'Fontsize', 12)
    grid on;
    box on;
    xlabel('Current density, {\it{j}} [A/cm^2]');
    ylabel('cell voltage, {\it{v}} [V]');
    
  end
  
end
%%                                    FIT POLARIZATION CURVES

% INPUT DATA
A_cell    = 238;            % cell active area       [cm2]
N_cell    = 178;            % number of cells    
p_atm     = 1;              % atmosferic pressure    [atm]
p_an      = p_atm * 3.50;   % anode pressure         [atm]
p_cat     = p_atm * 3.50;   % cathode pressure       [atm] 
rH_H2     = 0.54;           % H2 relative humidity
rH_O2     = 0.54;           % O2 relative humidity
z_m       = 30e-04;         % membrane dry thickness [cm]
lambda_H2 = 1.3;
lambda_O2 = 1.5;
F         = 96485;          % Faraday's constant     [C / mol]
R         = 8.314;          % Universal gas constant [J / mol-K]
    
% FITTING
[xi, fval] = fmincon(@(xi) PEMFC_H2O2polfitting(I, V, T, A_cell, N_cell, F, lambda_H2, lambda_O2, ...
  rH_H2, rH_O2, p_an, p_cat, z_m, xi), [-0.949, 0.003540, 7.6e-05, -1.93e-04, 3.3e-03, 0.634, 1.81, 2]);
xi1 = xi(1);
xi2 = xi(2);
xi3 = xi(3);
xi4 = xi(4);
xi5 = xi(5);
xi6 = xi(6);
xi7 = xi(7);
xi8 = xi(8);
disp('---------------------------------------------------------------------------------------------------------')
disp('                          Approximating coefficients - Linear')
disp('---------------------------------------------------------------------------------------------------------')
Table = table(xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,'RowNames', {'xi'});
disp(Table);
disp('---------------------------------------------------------------------------------------------------------')
fprintf('MSE = %4.4f \n', fval)
disp('---------------------------------------------------------------------------------------------------------')

% PLOT RESULTS
for i = 1 : length(T)

  a1 = fieldnames(I);
  a2 = char(a1(i));
  j  = extractfield(I, a2)';
  b1 = fieldnames(V);
  b2 = char(b1(i));
  U  = extractfield(V, b2)';  
    
  % GAS RELATIVE HUMIDITY
  % Water saturation pressure [atm]
  pH2Osat = WaterSaturationPressure (T(i));

  % Water molar fraction in fed gas
  xH2O_H2 = rH_H2 * pH2Osat / p_an;
  xH2O_O2 = rH_O2 * pH2Osat / p_cat;

  % Current density [A / cm2]
  J = j * A_cell;

  % Inlet molar flow [mol/s]
  nH2in     = J * A_cell * N_cell / (2 * F) * lambda_H2;
  nO2in     = J * A_cell * N_cell / (4 * F) * lambda_O2;
  nH2Oin_H2 = nH2in * xH2O_H2 / (1 - xH2O_H2);
  nH2Oin_O2 = nO2in * xH2O_O2 / (1 - xH2O_O2);

  % Reacting molar flow [mol/s]
  nH2r     = J * A_cell * N_cell / (2 * F);
  nO2r     = J * A_cell * N_cell / (4 * F);
  nH2Or_H2 = J * A_cell * N_cell / (2 * F);

  % Outlet molar flow [mol/s]
  nH2out     = nH2in - nH2r;
  nO2out     = nO2in - nO2r;
  nH2Oout_H2 = nH2Oin_H2 + nH2Or_H2;
  nH2Oout_O2 = nH2Oin_O2;

  % Inlet partial pressure [atm]
  p_H20  = (1 - xH2O_H2) * p_an;
  p_O20  = (1 - xH2O_O2) * p_cat;

  % Outlet partial pressure [atm]
  pH2out     = nH2out ./ (nH2out + nH2Oout_H2) * p_an;
  pO2out     = nO2out ./ (nO2out + nH2Oout_O2) * p_cat;

  % Average partial pressure [atm]
  pH2avg = (p_H20 + pH2out) * 0.35;
  pO2avg = (p_O20 + pO2out) * 0.35;

  % Open circuit voltage [V]
  E0 = 1.02 - 0.85e-03 * (T(i) - 298.15) + 4.3085e-05 * T(i) * (log(pH2avg) + log(pO2avg) / 2); 

  % Activation losses [V]
  c_O2   = pO2avg / (5.08e06 * exp(-498 / T(i)));  
  V_act0 = - (xi(1) + xi(2) * T(i) + xi(3) * T(i) * log(c_O2));
  V_act  = (V_act0 - xi(4) * T(i) * log(J));

  % Ohmic losses [V]
  psi   = 23;
  phi1  = 1 + 0.03 * j + 0.062 * (T(i) / 333)^2 * j.^2.5;
  phi2  = psi - xi(6) - 3 * j;
  phi3  = 4.18 * (T(i) - 333) / T(i); 
  rho_m = xi(7) * 100 * phi1 ./ (phi2 .* exp(phi3));   
  R_m   = rho_m * z_m / A_cell;     
  R_eq  = xi(5) + R_m;
  V_ohm = J .* R_eq;

  % Diffusion losses [V]
  B     = 0.0050;          
  j_max = xi(8);     
  V_dif = - B * log(1 - j / j_max);

  % Voltage [V]
  u = E0 - V_dif - V_act - V_ohm;   
    
  % ERROR CALCULATION
  res = sum(1/length(u) * (U - u).^2);
  
  figure;
  plot(j, U, j, u, 'Linewidth', 2)
  set(gca, 'Fontsize', 12)
  grid on;
  box on;
  xlabel('Current density, {\it{j}} [A/cm^2]');
  ylabel('cell voltage, {\it{v}} [V]');
    
end