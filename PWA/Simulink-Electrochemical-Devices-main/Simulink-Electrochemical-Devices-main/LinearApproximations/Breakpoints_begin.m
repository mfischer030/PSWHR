clear variables;
close all;
clc;
D   = load('SOFC_data.txt');
q   = D(:, 1)';
p   = D(:, 4)';
r   = D(:, 5)';
Nbp = 1;

plot(q, p./q)

breakpoints = optimal_pwa_fit(Nbp, q, p, 0, 0);
fc_x_breakpoint_vals = breakpoints(1:end/2);
fc_y_breakpoint_vals = breakpoints(end/2+1:end);
[~, fc_pwa_slopes, fitted_y] = breakpoint_eval_fun([fc_x_breakpoint_vals; fc_y_breakpoint_vals], q, p, 0);