%--------------------------------------------------------------------------------------------------%
%                                             Heat pump                                            %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, July 2016

%--------------------------------------------------------------------------------------------------%
%                                    Fitting of size parameters                                    %
%--------------------------------------------------------------------------------------------------%

function res = fitting_function(IN, OUT, S, a, k)

% INITIALIZATION
res = 0;

% LINEAR FITTING
for i = 1 : length(S) 
  in  = (k(4) * S(i) + k(5)) .* (1 - a{i}) + (k(6) * S(i) + k(7)) .* a{i};
  out = k(1) * in + k(2) * S(i) + k(3);

  res1 = 1/length(IN{i}) * ((IN{i} - in).^2 / S(i));
  res2 = 1/length(OUT{i}) * ((OUT{i} - out).^2 / S(i));
  res  = res + sum(res1 + res2);
end

end