%-------------------------------------------------------------------------%
%                    PEMFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                             Process Engineering Institute, September 2015


%-------------------------------------------------------------------------%
%          Grey-box Steady-State PEMFC Design - Polarization Curve        %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Preliminary design of a proton exchange membrane fuel cell (PEMFC). The 
% cell is modeled as a grey box component, where physical laws are combined
% with data from literature. Static and dynamic behavior are described, and
% both electrochemical and thermal features of the PEM fuel cell are
% captured.
% Modeling parameter were chosen from the literature, by assuming the 
% following:
% 1. 0-D model: greybox approach;
% 3. Ideal gases;
% 4. Dynamic describing a H2-fuelled channel;
% 5. Negligible pressure drop along the channels;
% 6. Constant operating temperature (design conditions).

%%
function [V, P] = PEMFC_cellstack(I, N_cell, A, T, p_H2, p_O2, p_H2O)

%%                      ELECTROCHEMICAL MODEL      

% CONSTANTS
F = 96485;             % Faraday's constant                                [C / mol]
R = 8.314;             % Universal gas constant                            [J / mol-K]

% DESIGN CURRENT DENSITY [A / cm2]
i = I / A;

% NERNST MODIFIED POTENTIAL [V]
E00 = 1.229 + (R * T) / (2 * F) * log(p_H2 * p_O2^0.5 / p_H2O);
E0  = 1.229 - 0.85e-03 * (T - 298.15) + 4.3085e-05 * T * ...
 (log(p_H2) + log(p_O2) / 2);      % Open-circuit potential   

E   = E0;                          % Open-circuit (reversible) potential 

% ACTIVATION LOSSES
% Concentration at membrane interface - Data from [1]
c_H2 = p_H2 / (2.50e05 * exp(-498 / T));   % H2 concentration at TPB       [mol / cm3]
c_O2 = p_O2 / (5.08e06 * exp(-498 / T));   % O2 concentration at TPB       [mol / cm3]

% Fitting parameters - Data from [1]
%xi2 = 0.00286 + 0.0002 * log(A) + 4.3e-05 * log(c_H2);
xi1 = -0.949;
xi2 = 0.003420;
xi3 = 7.6e-05;
xi4 = -1.93e-04;

% Overvoltage [V]
V_act0 = - (xi1 + xi2 * T + xi3 * T * log(c_O2));
V_act  = (V_act0 - xi4 * T * log(I)) / 2.85;

% OHMIC LOSSES
% Fitting parameters - Data from [1]
xi5 = 0.0001;   % Equivalent contact resistance to electrons conduction   [ohm]
psi  = 23;
phi1 = 1 + 0.03 * i + 0.062 * (T / 333)^2 * i^2.5;
phi2 = psi - 0.634 - 3 * i;
phi3 = 4.18 * (T - 333) / T;

% Membrane resistance to proton conduction  
rho_m = 181.6 * phi1 / (phi2 * exp(phi3));   % Membrane resistivity        [ohm cm]
z_m   = 178e-04;                             % Membrane thickness          [cm]
R_m   = rho_m * z_m / A;     % Equiv. membrane resistance to proton cond.  [ohm]

% Equivalent internal resistance [ohm]
R_eq = xi5 + R_m;       

% Overvoltage [V]
V_ohm = I * R_eq;

% CONCENTRATION POLARIZATION LOSSES
% Fitting parameters - Data from [2]
B     = 0.012;          
i_max = 1.5;     % Limit current density                                   [A / cm2]     
% B     = 0.016;          
% i_max = 1.2;     % Limit current density                                   [A / cm2]     

% Overvoltage [V]
V_conc = - B * log(1 - i / i_max);

% DESIGN VOLTAGE
V_ = E - V_conc - V_act - V_ohm;   % Single cell voltage                   [V]
V  = V_ * N_cell;                  % Stack voltage                         [V]

% DESIGN POWER
P      = V * I;                      % Stack electric power                [W]

end