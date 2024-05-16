%-------------------------------------------------------------------------%
%                     SOFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                               Process Engineering Institute, October 2015


%-------------------------------------------------------------------------%
%           Grey-box Steady-State SOFC Design - Polarization Curve        %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Preliminary design of a solid oxide fuel cell (SOFC). The cell is modeled 
% as a 0D grey box component, where physical laws are combined with data 
% from literature. Static and dynamic behavior are described. Only 
% electrochemical features of the PEM fuel cell are captured.
% Modeling parameter were chosen from the literature, by assuming the 
% following:
% 1. 0-D model: greybox approach;
% 2. Ideal gases;
% 3. Dynamic of the stack describing a H2-fuelled channel;
% 4. Negligible pressure drop along the channels;
% 5. Constant operating temperature (design conditions);
% 6. Complete reforming of methane and shift reaction at equlibrium
% 7. CO not electrochemically oxidized, but converted through the shift
%    reaction

%%
function [V_, P, eta_e] = SOFC_stack(I, p_atm, p_an, p_cat, T, yy, xx, ...
    N_cell, A_cell, A_a, A_c, A_el, B_a, B_c, B_el, alpha, E_act_a, ...
    E_act_c, i_max, gamma_a, gamma_c, t_a, t_c, t_el, theta, F, R, U, ...
    lambda, LHV_CH4, LHV_CO, LHV_H2, M_CH4, M_CO, M_H2)

%%                      ELECTROCHEMICAL MODEL      

% PRE-ALLOCATION
nn_ain = zeros(1, length(yy));
nn_cin = zeros(1, length(xx));
nn_aav = zeros(1, length(yy));
nn_cav = zeros(1, length(xx));
yy_av  = zeros(1, length(yy));
xx_av  = zeros(1, length(xx));
pp_aav = zeros(1, length(yy));
pp_cav = zeros(1, length(xx));

% CURRENT DENSITY [A / cm2]
i = I / (A_cell);

% CHEMICAL REACTIONS AND MOLAR FLOW RATES 
%--------------------------------------------------------------------------
% The reaction occurring in the fuel cell are the following:
%  Methane reforming:
%   1. CH4 + H20 ---> CO + 3H2              (xi1)
%   2. CO + H20 ---> CO2 + H2               (xi2)

%  Electrochemical reactions:
%   Anode:     H2 + O2mm ---> H20 + 2e
%   Cathode:   1/2 O2 + 2e ---> O2mm 
%  
%   3. H2 + 1/2 O2 ---> H2O                 (xi3)
%--------------------------------------------------------------------------

% Inlet flow rates [mol / s]
n_a = I * N_cell / (2 * F * U * (4 * yy(1) + yy(2) + yy(6)));  
n_c = I * N_cell * lambda / (4 * F * xx(5));                

% Reaction rates, xi, [mol / s]
xi1 = n_a * yy(1);                                   % Reaction 1             
xi3 = U * n_a * (4 * yy(1) + yy(2) + yy(6));         % Reaction 3     
Ks  = 0.0126 * exp (4639 / T);                       % Eq. cons. from [5]    
res = @(xi2) Ks - ((n_a * yy(3) + xi2) * ...
      (n_a * yy(6) + 3 * xi1 + xi2 - xi3))...
      / ((n_a * yy(2) + xi1 - xi2) * (n_a * yy(7) - xi1 - xi2 + xi3));
opt = optimset('FunValCheck', 'on', 'Display', 'off');
xi  = fsolve(res, n_a * (yy(1) + yy(2)) / 2, opt);
xi2 = xi;                                            % Reaction 2

% Inlet flow rates - Anode side [mol / s]
for j = 1 : length(yy)
 nn_ain(j) = n_a * yy(j);
end

% Outlet flow rates - Anode side [mol / s]
nn_aout(1) = n_a * yy(1) - xi1;                 
nn_aout(2) = n_a * yy(2) + xi1 - xi2;  
nn_aout(3) = n_a * yy(3) + xi2;                        
nn_aout(4) = n_a * yy(4);                     
nn_aout(5) = n_a * yy(5); 
nn_aout(6) = n_a * yy(6) + 3 * xi1 + xi2 - xi3;  
nn_aout(7) = n_a * yy(7) - xi2 - xi1 + xi3; 
nn_aout(8) = n_a * yy(8);                        

% Inlet flow rates - Cathode side [mol / s]
for j = 1 : length(xx)
 nn_cin(j) = n_c * xx(j);
end  

% Outlet flow rates - Cathode side [mol / s]
nn_cout(1) = n_c * xx(1); 
nn_cout(2) = n_c * xx(2);  
nn_cout(3) = n_c * xx(3);
nn_cout(4) = n_c * xx(4); 
nn_cout(5) = n_c * xx(5) - 0.5 * xi3;                                            

% Average molar flow rates - Anode side [mol / s]
for j = 1 : length(yy)
 nn_aav(j) = (nn_ain(j) + nn_aout(j)) / 2;
end

% Average molar flow rates - Cathode side [mol / s]
for j = 1 : length(xx)
 nn_cav(j) = (nn_cin(j) + nn_cout(j)) / 2;
end

% Average molar fractions - Anode side [mol / s]
for j = 1 : length(yy)
 yy_av(j) = nn_aav(j) / sum(nn_aav);
end

% Average molar fractions - Cathode side [mol / s]
for j = 1 : length(xx)
 xx_av(j) = nn_cav(j) / sum(nn_cav);
end

% NERNST VOLTAGE
% Partial pressures (Dalton's law) - Anode side [Pa]
for j = 1 : length(yy)
 pp_aav(j) = yy_av(j) * p_an;
end

% Partial pressures (Dalton's law) - Cathode side [Pa]
for j = 1 : length(xx)
 pp_cav(j) = xx_av(j) * p_cat;
end

% Ideal Reversible voltage [V]
E0 = 1.2723 - 2.7645e-4 * T;

% Open-circuit voltage [V]
dummyA = (R * T) / (2 * F);
E      = E0 + dummyA * log ((pp_aav(6) * (pp_cav(5)^(0.5))) / pp_aav(7));
E      = E * theta; 

% ACTIVATION VOLTAGE LOSSES
% Activation current density [A / cm2]
i0_a = gamma_a * ((pp_aav(6) / p_atm) * ...
    (pp_aav(6) / p_atm)) * ...
    (exp (- E_act_a / (R * T)));                        % Anode side
i0_c = gamma_c * ((pp_cav(5) / p_atm)^(0.25)) * ...
    (exp (- E_act_c / (R * T)));                        % Cathode side
    
% Activation parameters
K1 = R * T / F;
K2 = 1 / (2 * i0_a);
K3 = 1 / (2 * i0_c);

% Voltage losses
V_act_a = K1 * asinh(K2 * i) * alpha / 0.5;    % Anode activation losses   [V]
V_act_c = K1 * asinh(K3 * i) * alpha / 0.5;    % Cathode activation losses [V]            
V_act   = V_act_a + V_act_c;                   % Total activation losses   [V]        

% OHMIC LOSSES
% Ohmic resistivity
rho_a  = A_a * exp (B_a / T);     % Anode resistivity                      [ohm cm]       
rho_c  = A_c * exp (B_c / T);     % Cathode resistivity                    [ohm cm]    
rho_el = A_el * exp (B_el / T);   % Electrolyte resistivity                [ohm cm]   

% Overall resisitivity [ohm cm2]
rho = rho_a * t_a + rho_c * t_c + rho_el * t_el; 

% Voltage losses [V]
V_ohm = i * rho;

% CONCENTRATION VOLTAGE LOSSES [V]
V_con = abs(dummyA * log(1 - i / i_max));

% VOLTAGE 
V_ = E - V_act - V_ohm - V_con;    % Cell voltage                          [V]
V  = V_ * N_cell;                  % Stack voltage                         [V]

% ELECTRICAL POWER OUTPUT [W]
P = V * I;

% THERMAL POWER INPUT [W]
Th = n_a * (LHV_CH4 * yy(1) * M_CH4 + LHV_CO * yy(2) * M_CO + ...
    LHV_H2 * yy(6) * M_H2);

% ELECTRICAL EFFICIENCY
eta_e = P / Th;

end
