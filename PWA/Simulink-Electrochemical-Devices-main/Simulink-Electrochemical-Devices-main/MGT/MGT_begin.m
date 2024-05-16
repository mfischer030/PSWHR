%--------------------------------------------------------------------------------------------------%
%                                               IMES                                               %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, June 2016

%--------------------------------------------------------------------------------------------------%
%                                               MGT                                                %
%--------------------------------------------------------------------------------------------------%

%%                                 DESCRIPTION AND ASSUMPTIONS

% The thermodynamic model is based on the Capstone C30 micro-gas turbine, a 30 kW electrical power
% output MGT. The model is based on data from literature along with some general assumptions.

% The turbine performance are build by finding possible working points for which the operating 
% conditions of the compressor match the operating conditions of the expander. Typical 
% characteristic maps for a MGT are used for both the components to define the working points. The
% points giving the maximum efficiency are selected and define the conversion performance.

% A few assumptions are added:
% 1. Single shaft MGT: same rotational speed for both compressor and turbine
% 2. Minimum threshold for expander power for being able to carry the compressor
% 3. --

% The model was originally developed by Ilias Papagiannis (LEC, D-MAVT, ETH)

%%                                       INITIALIZATION

clear all
close all
clc;

%%                                         INPUT DATA

load('Capstone_Oper_opt.mat')
load('Capstone_data.mat')

% ANALYSIS
flag_plot   = 1;    % 0: no plot; 1: with plot
flag_appr   = 0;    % 0: no approximation, 1: PWA approximations of conversion curves
approx_size = 1;

%%                                     FITTING EFFICIENCY

% NORMALIZED INLET POWER
P_MGT  = linspace(0.08,1,50);
C_MGT    = 31.4;

% ELECTRICAL EFFICIENCY
fit1 = Oper_opt(:,1) / 30;
fit2 = Oper_opt(:,2);
[funct_eta, gof_Eff] = MGT_exp_fit(fit1, fit2);

etae_MGT = funct_eta(P_MGT)';

% ELECTRICAL POWER
P_MGTin = P_MGT ./ etae_MGT;

% TURBINE EXIT TEMPERATURE
fit3 = Oper_opt(:,1);
fit4 = Oper_opt(:,5);
[xData, yData] = prepareCurveData(fit3, fit4);
fit_TET = fittype( 'poly1' );
[funct_TET, gof_TET] = fit(xData, yData, fit_TET);

TET = funct_TET(P_MGT);

% TURBINE EXIT FLOW
fit5 = Oper_opt(:,1);
fit6 = Oper_opt(:,7);
[xData, yData] = prepareCurveData(fit5, fit6);
fit_mair = fittype( 'poly1' );
[funct_mexh, gof_mair] = fit(xData, yData, fit_mair);

mair = funct_mexh(P_MGT);

% PRODUCED THERMAL POWER
Thw      = 45 + 273.15;    % Temperature of inlet hot water [K]
DeltaThw = 5;              % Heat exchanger DeltaT          [K]

Q_MGT     = mair * 1.003 .* (TET - Thw - DeltaThw);
etat_MGT  = Q_MGT' ./ (P_MGTin * C_MGT);
q_MGT     = etat_MGT .* P_MGTin;

% ELECTRICITY-TO-HEAT RATIO
% R_MGT = etae_MGT(end) / etat_MGT(end);
R_MGT = P_MGT ./ (etat_MGT .* P_MGTin);

%%                                   PLOT RESULTS
if flag_plot == 1;
  figure;
  subplot(1,3,1);
  plot(P_MGTin, P_MGT, 'Linewidth', 1.5);
  subplot(1,3,2);
  plot(P_MGTin, q_MGT, 'Linewidth', 1.5);
  subplot(1,3,3);
  plot(P_MGT, q_MGT, 'Linewidth', 1.5);
end

% REDUCED ORDER MODELS
if flag_appr == 1
  
  % Time constants [hr]
  tau_P = 10 / 3600;
  tau_Q = 50 / 3600;
  
  % PWA coefficients 1 kW
  N_bp     = 4;
  idx      = round(linspace(1, 50, N_bp));
  x_bp_val = (P_MGTin(idx)' - 0.6010);
  y_bp_val = (P_MGT(idx)');
  figure;
  plot(P_MGTin, P_MGT,'o', x_bp_val, y_bp_val,  'Linewidth', 1.5);
  m_etae   = zeros(1,length(x_bp_val)-1);
  q_etae   = zeros(1,length(x_bp_val)-1);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,1) = m_etae;
  qq(:,1) = q_etae;
  
  % Power to heat ratio
  r = polyfit([P_MGT(1) P_MGT(end)], [q_MGT(1), q_MGT(end)], 1);
  
  % PWA coefficients 100 kW
  P_MGT_   = P_MGT * 100;
  P_MGTin_ = (P_MGTin - 0.6010) * 100;
  x_bp_val = P_MGTin_(idx)';
  y_bp_val = P_MGT_(idx)';
  plot(P_MGTin_,P_MGT_,'o', x_bp_val, y_bp_val, 'Linewidth', 1.5);
  for i = 1 : length(x_bp_val) - 1
    m_etae(i) = (y_bp_val(i+1) - y_bp_val(i)) / (x_bp_val(i+1) - x_bp_val(i));
    q_etae(i) = y_bp_val(i+1) - m_etae(i) * x_bp_val(i+1);
  end
  mm(:,2) = m_etae;
  qq(:,2) = q_etae;
  
  polm = zeros(size(mm,1),1);
  polq = zeros(size(mm,1),2);
  
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit([1 100], mm(i,:), 0);
    polq(i,:) = polyfit([1 100], qq(i,:), 1);
  end

  % SAVING RESULTS
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech');
  Tech.MGT.etal  = 0.260;
  Tech.MGT.alpha = polm;
  Tech.MGT.beta  = polq;
  Tech.MGT.taue  = tau_P;
  Tech.MGT.taut  = tau_Q;
  Tech.MGT.R     = r;
  Tech.MGT.aux   = [0.6010, 0.08];
  save('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Tech', 'Tech');
end

eta_MGT     = [0.85 1.00, 1.11, 1.26, 1.26, 1.26];
C_MGT       = [1, 30, 65, 200, 600, 800];
N_bp        = 4;
idx         = round(linspace(1, 50, N_bp));

if approx_size == 1
  
  for i = 1 : length(eta_MGT)
    
    m_etae      = zeros(1,N_bp-1);
    q_etae      = zeros(1,N_bp-1);
    mm          = zeros(N_bp-1, length(C_MGT));
    qq          = zeros(N_bp-1, length(C_MGT));

    % NORMALIZED INLET POWER
    P_MGT  = linspace(0.08,1,50) * C_MGT(i);

    % ELECTRICAL POWER
    P_MGTin = P_MGT ./ (etae_MGT * eta_MGT(i));
    
    % PWA coefficients
    m1       = (P_MGT(17) - P_MGT(1)) / (P_MGTin(17) - P_MGTin(1));
    q1       = P_MGT(17) - m1 * P_MGTin(17);
    x1(i)    = -q1 / m1;
    x_bp_val = (P_MGTin(idx)' - x1(i));
    y_bp_val = P_MGT(idx)';
    plot(P_MGTin,P_MGT,'o', x_bp_val, y_bp_val, 'Linewidth', 1.5);
    for j = 1 : length(x_bp_val) - 1
      m_etae(j) = (y_bp_val(j+1) - y_bp_val(j)) / (x_bp_val(j+1) - x_bp_val(j));
      q_etae(j) = y_bp_val(j+1) - m_etae(j) * x_bp_val(j+1);
    end
    mm(:,i) = m_etae;
    qq(:,i) = q_etae;
  
  end
  
  % Coefficient calculation
  polm = zeros(size(mm,1),1);
  polq = zeros(size(mm,1),2);
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit(C_MGT, mm(i,:), 0);
    polm2(i,1) = mm(i,1) + 0.02;
    polq(i,:) = polyfit(C_MGT, qq(i,:), 1);
  end
  
end

if approx_size == 2
  
  % PARTIAL LOAD CONVERSION CURVES
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C30.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C65.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C200.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C600.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C800.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C1000.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C600dist.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C800dist.mat');
  load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C1000dist.mat');
  
  % COMPUTE VECTORS
  Q_C30   = eta_C30(:,1) ./ eta_C30(:,2);
  P_C30   = eta_C30(:,1);
  Q_C65   = eta_C65(:,1) ./ eta_C65(:,2);
  P_C65   = eta_C65(:,1);
  Q_C200  = eta_C200(:,1) ./ eta_C200(:,2);
  P_C200  = eta_C200(:,1);
  Q_C600  = eta_C600(:,1) ./ eta_C600(:,2);
  P_C600  = eta_C600(:,1);
  Q_C800  = eta_C800(:,1) ./ eta_C800(:,2);
  P_C800  = eta_C800(:,1);
  Q_C1000 = eta_C1000(:,1) ./ eta_C1000(:,2);
  P_C1000 = eta_C1000(:,1);
  Q.C30   = Q_C30;
  Q.C65   = Q_C65;
  Q.C200  = Q_C200;
  Q.C600  = Q_C600;
  Q.C800  = Q_C800;
  Q.C1000 = Q_C1000;
  P.C30   = P_C30;
  P.C65   = P_C65;
  P.C200  = P_C200;
  P.C600  = P_C600;
  P.C800  = P_C800;
  P.C1000 = P_C1000;
  S       = [30 65 200 600 800 1000];
  N_bp    = 2;
  mm      = zeros(N_bp-1, length(S));
  qq      = zeros(N_bp-1, length(S));
  pp      = zeros(1, length(S));

  for i = 1 : length(S)
    
    Q1 = fieldnames(Q);
    Q2 = char(Q1(i));
    P1 = fieldnames(P);
    P2 = char(P1(i));
  
    % INLET THERMAL POWER
    q  = extractfield(Q, Q2)'; 
    
    % OUTLET ELECTRICAL POWER
    p  = extractfield(P, P2)'; 
    
    % PWA COEFFICIENTS
    idx     = round(linspace(1, length(p), N_bp));
    x_bp_val = q(idx)';
    y_bp_val = p(idx)';
    for j = 1 : length(x_bp_val) - 1
      m_etae(j) = (y_bp_val(j+1) - y_bp_val(j)) / (x_bp_val(j+1) - x_bp_val(j));
      q_etae(j) = y_bp_val(j+1) - m_etae(j) * x_bp_val(j+1);
    end
    mm(:,i) = m_etae;
    qq(:,i) = q_etae;
    pp(i)   = q(1);
  
  end
  
  % Coefficient calculation
  polm = zeros(size(mm,1),1);
  polq = zeros(size(qq,1),2);
  for i = 1 : size(mm,1)
    polm(i,:) = polyfit(S, mm(i,:), 0);
    polq(i,:) = polyfit(S, qq(i,:), 1);
  end
  polp = polyfit(S, pp, 1);

  % PRINT COEFFICIENTS
  k1 = polm;
  k2 = polq(:,1);
  k3 = polq(:,2);
  k4 = [polp(:,1)];
  k5 = [polp(:,2)];
  k6 = [1];
  disp('--------------------------------------------------------------------------')
  disp('                   Approximating coefficients - PWA')
  disp('--------------------------------------------------------------------------')
  T = table(k1,k2,k3,k4,k5,k6,'RowNames', {'n=1'});
  disp(T);
  disp('--------------------------------------------------------------------------')
  
    for i = 1 : length(S)
    
      Q1 = fieldnames(Q);
      Q2 = char(Q1(i));
      P1 = fieldnames(P);
      P2 = char(P1(i));

      % INLET THERMAL POWER
      Q3  = extractfield(Q, Q2)'; 

      % OUTLET ELECTRICAL POWER
      P3  = extractfield(P, P2)'; 
      
      qmin = k4 * S(i) + k5;
      qmax = k6 * Q3(end);
      q    = linspace(qmin, qmax, length(Q3))'; 
      p    = k1 * q + k2 * S(i) + k3;
      R2   = 1/length(Q3) * sum((Q3 - q).^2 + (P3 - p).^2);
      figure;
      plot(Q3, P3, q, p, 'Linewidth', 2)
      set(gca, 'Fontsize', 12)
      grid on;
      box on;
      xlabel('Inlet thermal power, {\it{Q}} [kW]');
      ylabel('Outlet electrical power, {\it{P}} [kW]');
      title(['Conversion curve approximation for Capstone ', Q2])
      legend('Real efficiency curve', 'Approximate efficiency curve', 'Location', 'Northwest')
      annotation('textbox', [0.65, 0.125, 0.22, 0.08], 'String',['R^2 = ', num2str(R2)], ...
        'FontSize',10, 'Horizontalalignment', 'center', 'Verticalalignment', 'middle', ...
        'Background', 'white');
    end
    
end