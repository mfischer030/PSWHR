%-------------------------------------------------------------------------%
%                   IMES - Paolo Gabrielli - ETH Zurich                   %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                 Process Engineering Institute, March 2016

%-------------------------------------------------------------------------%
%            Function for breakpoints calculation in PWA fitting          %
%-------------------------------------------------------------------------%

function [ c_ineq, c_eq ] = breakpoint_slope_nlincon( breakpoints, input, output, floating_breakpoint, convexity )
% Determine whether slopes of line segments are monotone increasing or
% decreasing, depending on convexity.

c_eq = [];

x_vals = breakpoints(1:length(breakpoints)/2);

if (floating_breakpoint == 1),
    y_vals = breakpoints(length(breakpoints)/2+1:end);
else
    y_vals = interp1(input, output, x_vals);
end

slopes = [y_vals(1)/x_vals(1); diff(y_vals)./diff(x_vals)];
diff_vector = diff(slopes);

%  Should hull of interpolants b't breakpoints be convex or concave?
if (convexity == 0),
    c_ineq = diff_vector; 
else
    c_ineq = -diff_vector;
end


end

