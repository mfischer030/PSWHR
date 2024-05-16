%-------------------------------------------------------------------------%
%                   IMES - Paolo Gabrielli - ETH Zurich                   %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                 Process Engineering Institute, March 2016

%-------------------------------------------------------------------------%
%              Function for breakpoints optimal PWA fitting               %
%-------------------------------------------------------------------------%

function breakpoint = PWA_fit(num_breakpoints, input, output, ...
    floating_breakpoint)

% Compute optimal piecewise linear fit for relation between input and 
% output
initial_guess = [(input(1)+input(end))/(num_breakpoints+1) * ...
    (1:num_breakpoints) (output(1)+output(end))/(num_breakpoints+1) * ...
    (1:num_breakpoints)]';

% Each breakpoint must be larger than the previous
difference_matrix = -diag(ones(num_breakpoints-1,1),1) + ...
    eye(num_breakpoints);
difference_matrix = difference_matrix(1:end-1,:);
difference_matrix = [difference_matrix difference_matrix*0];
difference_limit  = zeros(num_breakpoints-1,1);

f = @(x)BreakpointsEval(x,input,output, floating_breakpoint);
options.MaxFunEvals = 9000;
options.Display     = 'off';
[breakpoint fval, exitflag, min_output] = fmincon(f,initial_guess, ...
    difference_matrix,difference_limit,[],[],[input(1) * ...
    ones(num_breakpoints,1); zeros(num_breakpoints,1)], ...
    [input(end)*ones(num_breakpoints,1);output(end) * ...
    ones(num_breakpoints,1)],[],options);

if (floating_breakpoint == 0),
 breakpoint(num_breakpoints+1:end) = interp1(input, output, ...
     breakpoint(1:num_breakpoints));
end

end