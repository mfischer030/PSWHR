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

function P_max = PV_IVcurve(k, T_ref, Ns, n_ref, q, Eg_ref, S_ref, alpha, beta, Voc_ref, ...
  Np, Ir_ref, I0_ref, Rs_ref, Rp_ref, Tair, S, flag_plot)

%%                                 FIVE-PARAMETER CALCULATION

% CELL TEMPERATURE
T = (Tair + (43 - 20) * S / 800) + 273.15;

% CELL IDEALITY FACTOR 
a_ref = k * T_ref * Ns * n_ref / q;
a     = a_ref * T / T_ref;

% BANDGAP ENERGY 
Eg = 1.17 - 0.000473 * (T^2 / (T + 636));

% CELL REVERSE SATURATION CURRENT
I0 = I0_ref * (T / T_ref)^3 * exp(q / (n_ref * k) * (Eg_ref / T_ref - Eg / T));

% CELL IRRADIANCE CURRENT
Ir = S / S_ref * (Ir_ref + alpha * (T - T_ref));

% CELL SHUNT RESISTANCE
Rp = S_ref / S * Rp_ref;

% CELL SERIES RESISTANCE
Rs = Rs_ref;

%%                                         I-V CURVE

% OPEN CIRCUIT VOLTAGE 
Voc = Voc_ref + beta * (T - T_ref);

% CURVES GENERATION
V   = linspace(0, Voc, 100)'; 
Idc = zeros(100, 1);
P   = zeros(100, 1);
options = optimset('Display','off');
for i = 1 : size(V,1)
  Idc(i) = fsolve(@(I) -I + Np * Ir - Np * I0 * (exp((V(i) + I * Ns / Np * Rs) / a) - 1) - ...
    (V(i) + I * Rs * Ns / Np) / (Ns / Np * Rp), 0, options);
  P(i) = Idc(i) * V(i);        
end
    
P_max = max(P);

if flag_plot == 1
  figure;
  subplot(1,2,1);
  plot(V, Idc, 'Linewidth', 1.5)
  subplot(1,2,2);
  plot(V, P, 'Linewidth', 1.5)
end

end