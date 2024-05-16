%% Micro Gas Turbine - Reduced Order Model
%% Ilias Papagiannis - 04/2016
%%

close all; clc; clear all

disp('---------------------------------------------')
disp('              Micro Gas Turbine              ')
disp('             Reduced Order Model             ')
disp('---------------------------------------------')
disp(' ')

%% Load the data of the model

% The results is a matrix containing all the operating points.
% Each column of the matrix contains the following quantities

% Data indices
% 1 = Net work
% 2 = Efficiency
% 3 = Rotational Spee
% 4 = TIT
% 5 = TET (Turbine Exhuast Temperature - Can be used for CHP)
% 6 = mass_fuel
% 7 = massflow (Massflow at the outlet of the machine - Used for CHP)

load('Capstone_Oper_opt.mat')
load('Capstone_data.mat')


%% Fit Efficiency

fit1=Oper_opt(:,1)/30;
fit2=Oper_opt(:,2);
[fitresult_Eff, gof_Eff] = MGT_exp_fit(fit1, fit2);
ft_Eff = fittype( 'exp2' );

% Plot fit with data.
figure(1);
h = plot(fitresult_Eff);
set(h,'Color','b','LineWidth',3)
hold on
textsize=16;
axis([-0.02 1.1 -0.005 0.35])
hx=xlabel('Normalized Electrical Power [-]');
set(hx,'FontSize',textsize,'FontWeight','Bold','FontName','arial')
hy=ylabel('Efficiency');
set(hy,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
titletext='Part Load Performance';
ht=title(char(titletext));
set(ht,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(gca,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
grid on
plot(Work_Capstone,Eff_Capstone/100,'r','Linewidth',3.0)
legend('Fitted curve','Capstone');

disp(' ')
disp('---------------------------------------------')
disp('         Exponential fit of Efficiency       ')
disp('---------------------------------------------')
disp(' ')
ft_Eff
fitresult_Eff



%% FIT TET
fit1=Oper_opt(:,1);
fit5=Oper_opt(:,5);

[xData, yData] = prepareCurveData( fit1, fit5 );

% Set up fittype and options.
ft_TET = fittype( 'poly1' );

% Fit model to data.
[fitresult_TET, gof_TET] = fit( xData, yData, ft_TET );

% Plot fit with data.
figure(2);
h=plot(fitresult_TET,xData,yData);
set(h,'Color','r','LineWidth',3)
hold on
plot(xData,yData,'.b','MarkerSize',30);
grid on

textsize=16;
axis([-0.5 33 475 565])
hx=xlabel('Net Power [kW]');
set(hx,'FontSize',textsize,'FontWeight','Bold','FontName','arial')
hy=ylabel('Turbine Exhaust Temperature [K]');
set(hy,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
titletext='Change of TET with power Output';
ht=title(char(titletext));
set(ht,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(gca,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(legend,...
    'Position',[0.2 0.703571428571432 0.232142857142857 0.120238095238095]);

disp(' ')
disp('---------------------------------------------')
disp('               Linear fit of TET             ')
disp('---------------------------------------------')
disp(' ')
ft_TET
fitresult_TET

%% FIT mass_fuel
fit1=Oper_opt(:,1);
fit6=Oper_opt(:,6)*1000;

[xData, yData] = prepareCurveData( fit1, fit6 );

% Set up fittype and options.
ft_mf = fittype( 'poly1' );

% Fit model to data.
[fitresult_mf, gof_mf] = fit( xData, yData, ft_mf );

% Plot fit with data.
figure(3);
h=plot(fitresult_mf,xData,yData);
set(h,'Color','r','LineWidth',3)
hold on
plot(xData,yData,'.b','MarkerSize',30);
grid on

textsize=16;
axis([-0.5 33 0.5 3.0])
hx=xlabel('Net Power [kW]');
set(hx,'FontSize',textsize,'FontWeight','Bold','FontName','arial')
hy=ylabel('Fuel consumption [g/s]');
set(hy,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
titletext='Change of m_f_u_e_l with power Output';
ht=title(char(titletext));
set(ht,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(gca,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(legend,...
    'Position',[0.2 0.703571428571432 0.232142857142857 0.120238095238095]);

disp(' ')
disp('---------------------------------------------')
disp('          Linear fit of Fuel Massflow        ')
disp('---------------------------------------------')
disp(' ')
ft_mf
fitresult_mf


%% FIT massflow
fit1=Oper_opt(:,1);
fit7=Oper_opt(:,7);

[xData, yData] = prepareCurveData( fit1, fit7 );

% Set up fittype and options.
ft_mair = fittype( 'poly1' );

% Fit model to data.
[fitresult_mair, gof_mair] = fit( xData, yData, ft_mair );

% Plot fit with data.
figure(4);
h=plot(fitresult_mair,xData,yData);
set(h,'Color','r','LineWidth',3)
hold on
plot(xData,yData,'.b','MarkerSize',30);
grid on

textsize=16;
axis([-0.5 33 0.10 0.35])
hx=xlabel('Net Power [kW]');
set(hx,'FontSize',textsize,'FontWeight','Bold','FontName','arial')
hy=ylabel('Air mass flow [kg/s]');
set(hy,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
titletext='Change of m_a_i_r with power Output';
ht=title(char(titletext));
set(ht,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(gca,'FontSize',textsize,'FontName','arial','FontWeight','Bold')
set(legend,...
    'Position',[0.2 0.703571428571432 0.232142857142857 0.120238095238095]);

disp(' ')
disp('---------------------------------------------')
disp('           Linear fit of Massflow            ')
disp('---------------------------------------------')
disp(' ')
ft_mair
fitresult_mair

clear fit1 fit2 fit5 fit6 fit7 h ht hx hy titletext textsize xData yData


