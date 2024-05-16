%--------------------------------------------------------------------------------------------------%
%                                       Micro Gas Turbine                                          %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, July 2016

%--------------------------------------------------------------------------------------------------%
%                            Size dependency of conversion efficiency                              %
%--------------------------------------------------------------------------------------------------%

clear variables
close all
clc

flag_plot_eff  = 0;
flag_plot_conv = 0;
flag_power = 'electrical';

%%                         PARTIAL LOAD EFFICIENCY CURVES

load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C30.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C65.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C200.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C600.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C800.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C1000.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C600dist.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C800dist.mat');
load('C:\Users\gapaolo\Documents\MATLAB\Energy Hub\Technologies\MGT\eta_C1000dist.mat');

%%                             PLOT EFFICIENCY CURVES

if flag_plot_eff  == 1
  
  figure;
  plot(eta.C30(:,1), eta.C30(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C30')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  figure;
  plot(eta.C65(:,1), eta.C65(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C65')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  figure;
  plot(eta.C200(:,1), eta.C200(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C200')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  figure;
  plot(eta.C600(:,1), eta.C600(:,2), 'Linewidth', 2);
  hold on;
  plot(eta_C600dist(:,1), eta_C600dist(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C600')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  legend('Maximum efficiency mode', 'Distributed efficiency mode', 'Location', 'Southeast')
  grid on;
  box on;
  hold off;

  figure;
  plot(eta.C800(:,1), eta.C800(:,2), 'Linewidth', 2);
  hold on;
  plot(eta_C800dist(:,1), eta_C800dist(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C800')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  legend('Maximum efficiency mode', 'Distributed efficiency mode', 'Location', 'Southeast')
  grid on;
  box on;
  hold off;

  figure;
  plot(eta.C1000(:,1), eta.C1000(:,2), 'Linewidth', 2);
  hold on;
  plot(eta_C1000dist(:,1), eta_C1000dist(:,2), 'Linewidth', 2);
  title('Partial load efficiency curve for Capstone C1000')
  xlabel('Outlet electrical power, {\it{P}} [kW]');
  ylabel('Electrical efficiency, {\it{\eta_e}} []');
  set(gca, 'Fontsize', 12);
  legend('Maximum efficiency mode', 'Distributed efficiency mode', 'Location', 'Southeast')
  grid on;
  box on;
  hold off;

end
%%                     PLOT INPUT-TO-OUTPUT CURVES FOR MAXIMUM EFFICIENCY MODE

% INPUT AND OUTPUT POWERS 
switch flag_power
  case 'electrical'
    index = 2;
  case 'thermal'
    index = 3;
end
Q_C30   = eta.C30(:,1) ./ eta.C30(:,2);
P_C30   = eta.C30(:,index) .* Q_C30;
Q_C65   = eta.C65(:,1) ./ eta.C65(:,2);
P_C65   = eta.C65(:,index) .* Q_C65;
Q_C200  = eta.C200(:,1) ./ eta.C200(:,2);
P_C200  = eta.C200(:,index) .* Q_C200;
Q_C600  = eta.C600(:,1) ./ eta.C600(:,2);
P_C600  = eta.C600(:,index) .* Q_C600;
Q_C800  = eta.C800(:,1) ./ eta.C800(:,2);
P_C800  = eta.C800(:,index) .* Q_C800;
Q_C1000 = eta.C1000(:,1) ./ eta.C1000(:,2);
P_C1000 = eta.C1000(:,index) .* Q_C1000;

if flag_plot_conv == 1
  
  % PLOT CONVERSION CURVES
  figure;
  plot(Q_C30, P_C30, 'Linewidth', 1);
  title('Conversion curve for Capstone C30')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  hold on;
  plot(Q_C65, P_C65, 'Linewidth', 2);
  title('Conversion curve for Capstone C30 - C65')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  hold on;
  plot(Q_C200, P_C200, 'Linewidth', 2);
  title('Conversion curve for Capstone C30 - C65 - C200')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  figure;
  plot(Q_C600, P_C600, 'Linewidth', 2);
  title('Conversion curve for Capstone C600')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  hold on;
  plot(Q_C800, P_C800, 'Linewidth', 2);
  title('Conversion curve for Capstone C30 - C800')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;

  hold on;
  plot(Q_C1000, P_C1000, 'Linewidth', 2);
  title('Conversion curve for Capstone C600 - C800 - C1000')
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  set(gca, 'Fontsize', 12);
  grid on;
  box on;
  
end

%%                              OPTIMAL REGRESSION OF EFFICIENCY APPROXIMATION

% LINEAR APPROXIMATION OF THE CONVERSION CURVES: P = k1 (k4 S a + k5 S (1 - a)) + k2 S + k3
%  P = outlet electrical power
%  k = fitting parameters
%  S = size
%  a = (Qmax - Q) / (Qmax - Qmin) = inlet thermal power

% COMPUTE VECTOR AND STRUCTS 
Q.C30   = Q_C30(round(linspace(1,length(Q_C30),29))');
Q.C65   = Q_C65(round(linspace(1,length(Q_C65),29))');
Q.C200  = Q_C200(round(linspace(1,length(Q_C200),29))');
Q.C600  = Q_C600(round(linspace(1,length(Q_C600),29))');
Q.C800  = Q_C800(round(linspace(1,length(Q_C800),29))');
Q.C1000 = Q_C1000(round(linspace(1,length(Q_C1000),29))');
P.C30   = P_C30(round(linspace(1,length(P_C30),29))');
P.C65   = P_C65(round(linspace(1,length(P_C65),29))');
P.C200  = P_C200(round(linspace(1,length(P_C200),29))');
P.C600  = P_C600(round(linspace(1,length(P_C600),29))');
P.C800  = P_C800(round(linspace(1,length(P_C800),29))');
P.C1000 = P_C1000(round(linspace(1,length(P_C1000),29))');
S       = [115.4 223.5 609.8 1829.3 2439.0 3048.8];
S1      = [115.4 223.5 609.8];
S2      = [609.8 1829.3 2439.0 3048.8];
a.C30   = 1 - (Q.C30(end) - Q.C30) ./ (Q.C30(end) - Q.C30(1));
a.C65   = 1 - (Q.C65(end) - Q.C65) ./ (Q.C65(end) - Q.C65(1));
a.C200  = 1 - (Q.C200(end) - Q.C200) ./ (Q.C200(end) - Q.C200(1));
a.C600  = 1 - (Q.C600(end) - Q.C600) ./ (Q.C600(end) - Q.C600(1));
a.C800  = 1 - (Q.C800(end) - Q.C800) ./ (Q.C800(end) - Q.C800(1));
a.C1000 = 1 - (Q.C1000(end) - Q.C1000) ./ (Q.C1000(end) - Q.C1000(1));

% FITTING
S   = S1;   % S1 for non-modular turbines: C30-C65-C200; S2 for modular turbines: C600-C800-C1000
idx = 0;    % 0 for non-modular turbines: C30-C65-C200; 2 for modular turbines: C600-C800-C1000

[k, fval] = fmincon(@(k) MGTsize_fitting(Q, P, S, a, idx, k), [0.33 -0.005 0 0.05 0 1 0 0 0 0]);
k1 = [k(1); k(8)];
k2 = [k(2); k(9)];
k3 = [k(3); k(10)];
k4 = [k(4); 0];
k5 = [k(5); 0];
k6 = [k(6); 0];
k7 = [k(7); 0];
disp('---------------------------------------------------------------------------------------------')
disp('                          Approximating coefficients - Linear')
disp('---------------------------------------------------------------------------------------------')
T = table(k1,k2,k3,k4,k5,k6,k7,'RowNames', {'n=1', 'n=2'});
disp(T);
disp('---------------------------------------------------------------------------------------------')
fprintf('MSE = %4.4f \n', fval)
disp('---------------------------------------------------------------------------------------------')

for i = 1 : length(S)
  j = i + idx;
  alpha  = fieldnames(a);
  alpha1 = char(alpha(j));
  alpha2 = extractfield(a, alpha1)';
  beta   = fieldnames(Q);
  beta1  = char(beta(j));
  beta2  = extractfield(Q, beta1)';  
  gamma  = fieldnames(P);
  gamma1 = char(gamma(j));
  gamma2 = extractfield(P, gamma1)';
  
  q      = (k(4) * S(i) + k(5)) * (1 - alpha2) + (k(6) * S(i) + k(7)) * alpha2;
  p      = k(1) * ((k(4) * S(i) + k(5)) .* (1 - alpha2) + (k(6) * S(i) + k(7)) .* alpha2) + k(2) * S(i) + k(3);
  res1 = 1/length(beta2) * ((beta2 - q).^2 / S(i));
  res2 = 1/length(gamma2) * ((gamma2 - p).^2 / S(i));
  R2  = sum(res1 + res2);
%   q = [27.8; q(4:end)];
%   p = [0.9; p(4:end)];
  figure;
  plot(beta2/gamma2(end), gamma2/gamma2(end), 'o', q/gamma2(end), p/gamma2(end), 'Linewidth', 2)
  set(gca, 'Fontsize', 12)
  grid on;
  box on;
  xlabel('Inlet thermal power, {\it{Q}} [kW]');
  ylabel('Outlet electrical power, {\it{P}} [kW]');
  title(['Conversion curve approximation for Capstone ', beta1])
  legend('Real efficiency curve', 'Approximate efficiency curve', 'Location', 'Northwest')
  annotation('textbox', [0.65, 0.125, 0.22, 0.08], 'String',['R^2 = ', num2str(R2)], ...
    'FontSize',10, 'Horizontalalignment', 'center', 'Verticalalignment', 'middle', ...
    'Background', 'white');
end
    
%%                                       HAKUNA MATATA
