%-------------------------------------------------------------------------%
%                     SOFC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                             Process Engineering Institute, September 2015


%-------------------------------------------------------------------------%
%                          SOFC-CHP Parameters                            %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic modeling of a solid oxide fuel cell (SOFC). 
% A lumped first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM fuel cell are captured.

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

%%                        THERMODYNAMIC MODELS

% STEADY-STATE POLARIZATION CURVE
if flag_dynam == 0      
 
 [V_des, P_des] = SOFC_cellstack (I_des, p_atm, p_an, p_cat, T, ...
     yy, xx, N_cell, A_cell, A_an, A_cat, A_el, B_an, B_cat, B_el, ...
     alpha, E_act_an, E_act_cat, i_max, gamma_an, gamma_cat, t_an, ...
     t_cat, t_el, theta, F, R, U, lambda, LHV_CH4, LHV_CO, LHV_H2, ...
     M_CH4, M_CO, M_H2);

% DESIGN CALCULATION
elseif flag_opcur == 0
 
 % Time horizon
 t  = 0 : 1 : (3600 * 3);
 z  = length(t);
 tf = t(z);
 
 nn = n_CH4in;
 % Inlet fuel time profile
 n_CH4in(1:3601)       = nn * 0.5;
 n_CH4in(3602:7201) = nn * 1.1;
 n_CH4in(7202:10801)   = nn * 0.10;
 
 % Simulation
 sim('SOFC_smlkCH4', t, [], [t; n_CH4in]');
 
 % Design efficiency calculation
 P    = PP.Data(end);
 Q    = QQ.Data(end);
 TT   = n_CH4in(1) * LHV_CH4 * M_CH4;
 etae = P / TT;
 etat = Q / TT;
 
% OFF-DESIGN CALCULATION
elseif flag_opcur == 1
 
 % Time horizon
 t  = 0 : 1 : 3600;
 z  = length(t);
 tf = t(z);
  
 % Pre-allocation for speeding
 nfuel = linspace(n_CH4in * 0.08, n_CH4in * 1.1, 50);
 w      = length(nfuel);
 P      = zeros(1, w);
 Q      = zeros(1, w);
 V      = zeros(1, w);
 J      = zeros(1, w);
 TT     = zeros(1, w);
 etae   = zeros(1, w);
 etat   = zeros(1, w);

 for k = 1 : w
   
  % Inlet fuel time profile
  n_CH4in(1:z) = nfuel(k);
   
  % Simulation
  run SOFC_CH4data
  sim('SOFC_smlkCH4', t, [], [t; n_CH4in]');
 
  % Design efficiency calculation
  TT(k)       = n_CH4in(1) * M_CH4 * LHV_CH4;
  P(k)        = PP.Data(end);
  Q(k)        = QQ.Data(end);
  V(k)        = VV.Data(end);
  J(k)        = II.Data(end);
  etae_aux(k) = P(k) / TT(k);
  etae(k)     = (P(k) + P_gamma) / TT(k);
  etat(k)     = Q(k) / TT(k);
 
 end
  
end