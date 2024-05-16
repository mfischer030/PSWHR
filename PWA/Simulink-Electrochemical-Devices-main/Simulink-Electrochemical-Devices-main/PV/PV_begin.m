%--------------------------------------------------------------------------------------------------%
%                                               IMES                                               %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, June 2016

%--------------------------------------------------------------------------------------------------%
%                                             PV panel                                             %
%--------------------------------------------------------------------------------------------------%

%%                                       INITIALIZATION

clear variables
close all
clc;

%%                                         INPUT DATA

% ANALYSIS
flag_plot = 0;                   % 0: no plot, 1: plot performance curve

% OPERATIONAL DATA - PS255 ECSOLAR Polycristallin
Voc_ref = 37.6;                  % reference open-circuit voltage                                   [V]
Isc_ref = 8.79;                  % reference short-circuit current                                  [A]
Vmp_ref = 30.4;                  % reference maximum power voltage                                  [V]
Imp_ref = 8.39;                  % reference maximum power current                                  [A]
S_ref   = 1000;                  % reference value of irradiance                                    [W/m2]
A       = 1640 * 992 * 1e-06;    % panel area                                                       [m2]

% OPERATIONAL DATA - BP365 Polycristallin ----------------------------------------------------------
% Voc_ref = 22.1;  
% Isc_ref = 3.99; 
% Vmp_ref = 17.6;  
% Imp_ref = 3.69;  
% S_ref   = 1000;   
% beta    = -0.08;        
% alpha   = 0.065*Isc_r;        
% --------------------------------------------------------------------------------------------------

% Temperature coefficients
Tnoct = 43;                      % [°C] 45+/-2
alpha = 0.055 * Isc_ref / 100;   % temperature coefficient of short-circuit voltage                 [1/K]
beta  = -0.33 * Voc_ref / 100;   % temperature coefficient of open-circuit voltage                  [1/K] 

% REFERENCE VALUES
Ns    = 60;                      % number of cells in series
Np    = 1;                       % number of cells in parallel
T_ref = 25 + 273.15;             % reference ambient temperature                                    [K]

% CONSTANT
q = 1.602e-19;                   % electronic charge                                                [C]
k = 1.3806503e-23;               % Boltzmann's constant                                             [J/K] 

%%                        STANDARD RATING CONDITIONS - TWO-PARAMETER MODEL

% PARAMETERS DEFINITION [2]
% Constant definition for avoiding cluttering
c1 = k * T_ref / q;
c2 = Voc_ref / (Ns * c1);
c3 = Vmp_ref / (Ns * c1);
c4 = Imp_ref / (Np * c1);
c5 = Isc_ref / (Np * c1);

% Exemplary +/- 10 cell temperature
T_ = 24 + 273.15;

% Bandgap energy for polycristalline silicon [eV]
Eg_ref = 1.17 - 0.000473 * ((T_ref)^2 / (T_ref + 636));
Eg     = 1.17 - 0.000473 * (T_^2 / (T_ + 636));

% Temperature-dependent factor for shunt conductance
KT = (T_ / T_ref)^3 * exp(Eg_ref * q / (k * T_ref) - Eg * q / (k * T_));


% SOLUTION OF THE NONLINEAR SYSTEM [2]
% Formulation of the domain of attraction: 0.5 <= n_ref <= 2, 0 <= Rs_ref <= Rs_max
Rs_max = (2 / c4 * (1 + lambertw(-exp(( c2 - 2 - 2 * c3) / 2))) + c3 / c4);

% Initial guess
x0 = [Rs_max/2 1];
lb = [0 0.5];
ub = [Rs_max 2];

% Two-equation set
foff = @(x) TwoParameters(x, Ns, k, T_ref, q, beta, T_, alpha, c1, c2, c3, c4, c5, KT);
options = optimset('TolFun', 1e-18, 'Display','iter');
y = lsqnonlin(foff,x0,lb,ub, options);

% Determination of reference cell resistance and cell ideality factor
Rs_ref = y(1);
n_ref  = y(2);

% Feasibility check
if Rs_ref == (c2 - c3) / c4;
  error('The cell resistance is not acceptable')
end

% ALGEBRAIC CALCULATION ON REMAINING PARAMETERS
% Reference irradiance current
Ir_ref = c1 * c2 * c4 / (c3 - c4 * Rs_ref) + ((c1 * c4 * (2 * c3 - c2)) / (c3 - c4 * Rs_ref) * ...
  (exp(c2 / n_ref) - 1 - c2 / n_ref * exp((c3 + c4 * Rs_ref) / n_ref))) / (exp(c2 / n_ref) + ...
  ((c4 * Rs_ref + c3 - c2) / n_ref - 1) * exp((c3 + c4 * Rs_ref) / n_ref));

% Reference cell reverse saturation current
I0_ref = c1 * c4 / (c3 - c4 * Rs_ref) * (2 * c3 - c2) / ((exp(c2 / n_ref) + ...
  ((c4 * Rs_ref + c3 - c2) / n_ref - 1) * exp((c3 + c4 * Rs_ref) / n_ref)));

% Reference shunt conductance
Gp_ref = (c4 * (((1 + (c3 - c4 * Rs_ref) / n_ref)) * exp((c3 + c4 * Rs_ref - c2) / n_ref) - 1)) ...
  / ((1 + ((c4 * Rs_ref + c3 - c2) / n_ref - 1) * exp((c3 + c4 * Rs_ref - c2) / n_ref)) * (c4 * ...
  Rs_ref - c3));

% Reference shunt resistance
Rp_ref = 1 / Gp_ref;

%%                                 I-V CURVE CALCULATION

% % INPUT DATA
% Tair = 25 + 273.15;   % ambient temperature                                                            [K]
% S = 1000;             % solar irradiance                                                               [W/m2]
% H = 8760;
% 
% % INITIALIZATION
% P_m = zeros(1, H);
% p_m = zeros(1, H);
% ETA = zeros(1, H);
% 
% % EFFICIENCY CALCULATION
% load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Data_Altstetten');
% 
% for i = 1 : size(Data.building.solarIrr,2)
%   for j = 1 : size(Data.building.solarIrr,1)
%     S    = Data.solarIrr(j,i) * 1000;
%     Tair = Data.Tair(j);
% 
%     if S > 10
%       % I-V curve and maximum power
%       P_max = PV_IVcurve(k, T_ref, Ns, n_ref, q, Eg_ref, S_ref, alpha, beta, Voc_ref, Np, ...
%         Ir_ref, I0_ref, Rs_ref, Rp_ref, Tair, S, flag_plot);
%       P_m(j,i) = P_max;
%       p_m(j,i) = P_max / A;
%       ETA(j,i) = p_m(j,i) / S;
%     else
%       p_m(j,i) = 0;
%       ETA(j,i) = 0;
%     end
%     fprintf('%4.2f%% percentage completed \n\n', j / (H * 1) * 100)
%   end
% end
% Data.PVefficiency = ETA;

%%                                                PLOT DEPENDENCIES

% INITIALIZATION
S = linspace (20,1000,50);
T = linspace(-10, 40, 50);
P_m = zeros(length(S), length(T));
p_m = zeros(length(S), length(T));
eta = zeros(length(S), length(T));

for i = 1 : length(S)
  s = S(i);
  
  for j = 1 : length(T)
    Tair = T(j);
    
    % I-V curve and maximum power
    P_max = PV_IVcurve(k, T_ref, Ns, n_ref, q, Eg_ref, S_ref, alpha, beta, Voc_ref, Np, ...
      Ir_ref, I0_ref, Rs_ref, Rp_ref, Tair, s, flag_plot);
    P_m(i,j) = P_max;
    p_m(i,j) = P_max / A;
    eta(i,j) = p_m(i,j) / s;

    fprintf('%4.2f%% percentage completed \n\n', (j + length(T) * (i - 1)) / length(S) / length(T) * 100)
    
  end
  
end
I = repelem(S',length(S));
T = repmat(T', length(T), 1);
e = reshape(eta', [length(S) * length(S), 1]);
etafit = fit([I, T], e, 'poly33');
plot(etafit, [I, T], e);
