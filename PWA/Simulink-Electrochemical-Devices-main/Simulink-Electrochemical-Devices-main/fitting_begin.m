%================================================================================================================================================
%                                               POWER AND SIZE DEPENDENCIES OF CONVERSION PERFORMANCE                                           %
%================================================================================================================================================

% Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
% Machanical Engineering, 
% Energy Science Center,
% Zürich, July 2016

%% DEFINE ANALYSIS

clear variables
close all
clc

% FLAGS
plot_input   = 0;              % 0) No plots       - 1) Plot efficiency curves
plot_results = 1;              % 0) No plots       - 1) Plot fitting results
PWA_points   = 0;              % 0) linear fitting - >0) number of breakpoints for piecewise affine (PWA) fitting

% TECHNOLOGY SELECTION
technology   = 'ALKEC';           % Fit the selected conversion technology
output       = 'thermal';   % Fit the selected output power: electrical, thermal, cooling
GT_type      = 'large';         % Fit the selected cathegory of gas turbines: aero, small, large

%% LOAD INPUT DATA

switch technology
  case 'GT'
    load([cd, '\GTperformance']);
    switch GT_type
      case 'aero'
        eta = GTperformance.aero;
      case 'small'
        eta = GTperformance.small;
      case 'large'
        eta = GTperformance.large;
    end
  case 'HP'
    load([cd, '\HPperformance']);
    eta = HPperformance;
end

%% PLOT EFFICIENCY CURVES

if plot_input == 1
  
  switch technology
    case {'GT'}
      for i = 1 : length(eta)
        figure;
        dummy = eta{i}.P(end) / 1000;
        subplot(2,1,1)
        plot(eta{i}.F, eta{i}.P);
        title(['Electrical efficiency curve GT - ', num2str(dummy), 'MW']);
        xlabel('Inlet thermal power, {\it{F}} [MW]');
        ylabel('Outlet electrical power, {\it{P}} [MW]');
        set(gca, 'Fontsize', 12);
        subplot(2,1,2)
        plot(eta{i}.F, eta{i}.Q);
        title(['Thermal efficiency curve GT - ', num2str(dummy), 'MW']);
        xlabel('Inlet thermal power, {\it{F}} [MW]');
        ylabel('Outlet thermal power, {\it{Q}} [MW]');
        set(gca, 'Fontsize', 12);  
      end
    case 'HP'
      for i = 1 : length(eta)
        figure;
        dummy = eta{i}.Qc(end);
        subplot(2,1,1)
        plot(eta{i}.Pc, eta{i}.Qc);
        title(['Partial load performance HP - ', num2str(dummy), 'kW - Cooling']);
        xlabel('Inlet electrical power, {\it{P}}_c [kW]');
        ylabel('Outlet cooling power, {\it{Q}}_c [kW]');
        set(gca, 'Fontsize', 12);
        subplot(2,1,2)
        plot(eta{i}.Ph, eta{i}.Qh);
        title(['Partial load performance HP - ', num2str(dummy), 'kW - Heating']);
        xlabel('Inlet electrical power, {\it{P}}_h [kW]');
        ylabel('Outlet heating power, {\it{Q}}_h [kW]');
        set(gca, 'Fontsize', 12);   
      end
  end
  
end
%% FITTING PROCEDURE

% INPUT AND OUTPUT POWERS 
S   = zeros(1, length(eta));
IN  = cell(1, length(eta));
OUT = cell(1, length(eta));
a   = cell(1, length(eta));

switch technology
  case {'GT'}
    varin = 'F';
    switch output
      case 'electrical'
        varout = 'P';
      case 'thermal'
        varout = 'Q';
    end
    for i = 1 : length(eta)
      S(i)   = eta{i}.F(end) / 1000;
      IN{i}  = eta{i}.(varin) / 1000;
      OUT{i} = eta{i}.(varout) / 1000;
      a{i}   = 1 - (IN{i}(end) - IN{i}) ./ (IN{i}(end) - IN{i}(1));
    end
  case 'HP'
    switch output
      case 'cooling'
        varin  = 'Pc';
        varout = 'Qc';
      case 'thermal'
        varin  = 'Ph';
        varout = 'Qh';
    end
    for i = 1 : length(eta)
      S(i)   = eta{i}.Ph(end);
      IN{i}  = eta{i}.(varin);
      OUT{i} = eta{i}.(varout) / 4 * 3.5;
      a{i}   = 1 - (IN{i}(end) - IN{i}) ./ (IN{i}(end) - IN{i}(1));
    end
end

% LINEAR APPROXIMATION OF THE CONVERSION CURVES: P = k1 (k4 S a + k5 S (1 - a)) + k2 S + k3
%  P = outlet electrical power
%  k = fitting parameters
%  S = size
%  a = (Fmax - F) / (Fmax - Fmin) = inlet thermal power

% FITTING FUNCTION
[k, fval] = fmincon(@(k) fitting_function(IN, OUT, S, a, k), [4.1 0 0 0.1 0 1 0]);
k1 = k(1);
k2 = k(2);
k3 = k(3);
k4 = k(4);
k5 = k(5);
k6 = k(6);
k7 = k(7);
disp('---------------------------------------------------------------------------------------------')
disp('                             Fitting coefficients - Linear')
disp('---------------------------------------------------------------------------------------------')
T = table(k1,k2,k3,k4,k5,k6,k7,'RowNames', {'k'});
disp(T);
disp('---------------------------------------------------------------------------------------------')
fprintf('MSE = %4.4f \n', fval)
disp('---------------------------------------------------------------------------------------------')
    
%% PLOT RESULTS

if plot_results == 1
  switch technology
    case {'GT'}
      for i = 1 : length(S)
        in  = (k(4) * S(i) + k(5)) .* (1 - a{i}) + (k(6) * S(i) + k(7)) .* a{i};
        out = k(1) * in + k(2) * S(i) + k(3);
        res1 = 1/length(IN{i}) * ((IN{i} - in).^2 / S(i));
        res2 = 1/length(OUT{i}) * ((OUT{i} - out).^2 / S(i));
        R2   = sum(res1 + res2);
        figure;
        plot(IN{i}, OUT{i}, in, out, 'Linewidth', 2)
        set(gca, 'Fontsize', 12)
        grid on;
        box on;
        xlabel('Inlet power, {\it{F}} [MW]');
        ylabel('Outlet power, {\it{P}} [MW]');
        title(['Conversion curve approximation - GT ', num2str(S(i)), ' MW'])
        legend('Real efficiency curve', 'Approximate efficiency curve', 'Location', 'Northwest')
        annotation('textbox', [0.65, 0.125, 0.22, 0.08], 'String',['MSE = ', num2str(R2)], 'FontSize',10, 'Horizontalalignment', 'center', ...
          'Verticalalignment', 'middle', 'Background', 'white');
      end
    case {'HP'}
      for i = 1 : length(S)
        in  = IN{i};
        out = k(1) * in + k(2) * S(i) + k(3);
        res1 = 1/length(IN{i}) * ((IN{i} - in).^2 / S(i));
        res2 = 1/length(OUT{i}) * ((OUT{i} - out).^2 / S(i));
        R2   = sum(res1 + res2);
        figure;
        plot(IN{i}, OUT{i}, in, out, 'Linewidth', 2)
        set(gca, 'Fontsize', 12)
        grid on;
        box on;
        xlabel('Inlet power, {\it{F}} [kW]');
        ylabel('Outlet power, {\it{P}} [kW]');
        title(['Conversion curve approximation - HP ', num2str(S(i)), ' kW'])
        legend('Real efficiency curve', 'Approximate efficiency curve', 'Location', 'Northwest')
        annotation('textbox', [0.65, 0.125, 0.22, 0.08], 'String',['MSE = ', num2str(R2)], 'FontSize',10, 'Horizontalalignment', 'center', ...
          'Verticalalignment', 'middle', 'Background', 'white');
      end
  end
end
