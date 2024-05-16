%-------------------------------------------------------------------------%
%                   IMES - Paolo Gabrielli - ETH Zurich                   %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                 Process Engineering Institute, March 2016

%-------------------------------------------------------------------------%
%            Function for breakpoints calculation in PWA fitting          %
%-------------------------------------------------------------------------%

function breakpoint = optimal_pwa_fit(num_breakpoints, input, output, floating_breakpoint, convexity)
% Compute optimal piecewise linear fit for relation between input and output

% Initial guess has breakpoints evenly spaced out
% initial_guess = [(input(1)+input(end))/(num_breakpoints+1)*(1:num_breakpoints) ...
%                  (output(1)+output(end))/(num_breakpoints+1)*(1:num_breakpoints)]';
scale = 0.6;
mean_x = (input(1)+input(end))/2;
range_x   = input(end)-input(1);
mean_y = (output(1)+output(end))/2;
range_y   = output(end) - output(1);
initial_guess = [((mean_x - scale*range_x/2):(scale*range_x/(num_breakpoints-1)):(mean_x + scale*range_x/2)) ...
                 ((mean_y - scale*range_y/2):(scale*range_y/(num_breakpoints-1)):(mean_y + scale*range_y/2))]';

% Each breakpoint must be larger than the previous
difference_matrix = -diag(ones(num_breakpoints-1,1),1)+eye(num_breakpoints);
difference_matrix = difference_matrix(1:end-1,:);
difference_matrix = [difference_matrix difference_matrix*0];
difference_limit  = zeros(num_breakpoints-1,1);

options.MaxFunEvals = num_breakpoints*1000;
options.Display = 'off';
f = @(x)breakpoint_eval_fun(x,input,output, floating_breakpoint);
g = @(x)breakpoint_slope_nlincon(x, input, output, floating_breakpoint, convexity);
[breakpoint, fval, exitflag, min_output] = fmincon(f,initial_guess,difference_matrix,difference_limit,[],[],...
                                                   [input(1)*ones(num_breakpoints,1); zeros(num_breakpoints,1)],...
                                                   [input(end)*ones(num_breakpoints,1);output(end)*ones(num_breakpoints,1)],...
                                                   g,options);

% If we specify a fixed breakpoint, then input y values have no effect
if (floating_breakpoint == 0),
    breakpoint(num_breakpoints+1:end) = interp1(input, output, breakpoint(1:num_breakpoints));
end

end