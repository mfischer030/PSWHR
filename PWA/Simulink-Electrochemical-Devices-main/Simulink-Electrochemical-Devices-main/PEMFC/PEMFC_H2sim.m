%-------------------------------------------------------------------------%
%                    PEMFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                             Process Engineering Institute, September 2015


%-------------------------------------------------------------------------%
%                         PEMFC-CHP Parameters                            %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic modeling of a proton exchange membrane fuel cell (PEMFC). 
% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM fuel cell are captured.

% The modeled based on a generic example of PEM fuel cell using methane 
% (1.4 kW). The polarization curve is based on data from PSI cell and from 
% literature.
% The operative power can range in 0-120% of the rated power.
% The fuel can be either natural gas or hydrogen. 
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

%%                        THERMODYNAMIC MODELS

% STEADY-STATE POLARIZATION CURVE
if flag_dynam == 0      
 
 [V_des, P_des] = PEMFC_cellstack(I_des, N_cell, A_cell, T, p_H2, p_O2, p_H2O);

% DESIGN CALCULATION
elseif flag_opcur == 0
 
 % Time horizon
 t  = 0 : 1 : 500;
 z  = length(t);
 tf = t(z);
 
 % Inlet fuel time profile
 n_H2in(1:z) = n_fuel;
 
 % Simulation
 sim('PEMFC_smlkH2', t, [], [t; n_H2in]');
 
 % Design efficiency calculation
 P    = PP.Data(end);
 Q    = QQ.Data(end);
 TT   = n_H2in(1) * M_H2 * LHV_H2;
 etae = P / TT;
 etat = Q / TT;
 
% OFF-DESIGN CALCULATION
elseif flag_opcur == 1
 
 % Time horizon
 t  = 0 : 1 : 500;
 z  = length(t);
 tf = t(z);
  
 % Pre-allocation for speeding
 nfuel  = linspace(n_fuel * 0.1, n_fuel * 1.25, 50);
 w      = length(nfuel);
 P      = zeros(1, w);
 Q      = zeros(1, w);
 TT     = zeros(1, w);
 etae   = zeros(1, w);
 etat   = zeros(1, w);

 for j = 1 : w
   
  % Inlet fuel time profile
  n_H2in(1:z) = nfuel(j);
   
  % Simulation
  run PEMFC_H2data
  sim('PEMFC_smlkH2', t, [], [t; n_H2in]');
 
  % Design efficiency calculation
  TT(j)       = n_H2in(1) * M_H2 * LHV_H2;
  P(j)        = PP.Data(end);
  Q(j)        = QQ.Data(end);
  etae_aux(j) = P(j) / TT(j);
  etae(j)     = (P(j) + P_gamma) / TT(j);
  etat(j)     = Q(j) / TT(j);
 
 end
  
end