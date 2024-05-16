%-------------------------------------------------------------------------%
%                   IMES - Paolo Gabrielli - ETH Zurich                   %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                 Process Engineering Institute, March 2016

%-------------------------------------------------------------------------%
%            Function for breakpoints calculation in PWA fitting          %
%-------------------------------------------------------------------------%

function [error_val, slope_data, fit_y_fin] = ...
    BreakpointsEval(breakpoints, input, output, floating_breakpoint)

plot_flag = 0;

if (nargin == 3),
 floating_breakpoint = 1; % breakpoints limited to linear interp b't pts
end

x_vals = breakpoints(1:length(breakpoints)/2);

if (floating_breakpoint == 1),
    
 y_vals = breakpoints(length(breakpoints)/2+1:end);
 
else
    
 y_vals = interp1(input, output, x_vals);

 if any(isnan(y_vals)),
     
  pleasePause = 1;
  
 end
 
end

error_val = 0;
slope_data = zeros(length(breakpoints)/2+1,1);

if (plot_flag == 1)
 hold off
 plot(input,output,'o-')
 hold on
 plot(x_vals,y_vals,'rx')
 grid on
end

fit_y1 = [];

for i = 1:length(x_vals),
 if (i == 1),
     
  % For first breakpoint, fit through origin
  x0 = 0;
  y0 = 0;
  x1 = x_vals(i);
  y1 = y_vals(i);

  cur_indices = find(input <= x_vals(i));
  
 else
     
  % For subsequent breakpoints, fit between adjacent breakpoints
  x0 = x_vals(i-1);
  y0 = y_vals(i-1);
  x1 = x_vals(i);
  y1 = y_vals(i);

  cur_indices = find((input >= x_vals(i-1)) & (input <= x_vals(i)));
  
 end

 line_slope = (y1-y0)/(x1-x0);
 slope_data(i) = line_slope;

 true_x = input(cur_indices);
 true_y = output(cur_indices);

 % Determine error from fitting
 fitted_y = y0 + line_slope*(true_x-x0);
 error_val = error_val + sum((fitted_y - true_y).^2);

 fit_y1 = [fit_y1 fitted_y];
 
 if (plot_flag == 1)
     
  plot(true_x,fitted_y,'ro')
  
 end


 if (i == length(x_vals))
     
  upper_indices = find(input > x_vals(end));

  true_x = input(upper_indices);
  true_y = output(upper_indices);

  x0 = x_vals(end);
  y0 = y_vals(end);
  line_slope = (true_x-x0)'\(true_y-y0)';
  slope_data(end) = line_slope;

  fitted_y = y0 + line_slope*(true_x-x0);
  error_val = error_val + sum((fitted_y - true_y).^2);

  if (plot_flag == 1)
      
    plot(true_x,fitted_y,'kx')
    
  end
  
 end

end

error_val = sqrt(error_val);

fit_y2 = fitted_y;
 
fit_y_fin = [fit_y1 fit_y2]; 

end