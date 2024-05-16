%--------------------------------------------------------------------------------------------------%
%                                       Micro Gas Turbine                                          %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, July 2016

%--------------------------------------------------------------------------------------------------%
%                                    Fitting of size parameters                                    %
%--------------------------------------------------------------------------------------------------%

function res = PEMFC_H2O2polfitting (I, V, T, A_cell, N_cell, F, lambda_H2, lambda_O2, ...
  rH_H2, rH_O2, p_an, p_cat, z_m, xi)

% SOLVING
res = 0;
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
  B     = 0.005;          
  j_max = xi(8);     
  V_dif = - B * log(1 - j / j_max);

  % Voltage [V]
  u = E0 - V_dif - V_act - V_ohm;   
    
  % ERROR CALCULATION
  res = res + sum(1/length(u) * (U - u).^2);
  
end

end