%--------------------------------------------------------------------------------------------------%
%                                       Micro Gas Turbine                                          %
%--------------------------------------------------------------------------------------------------%

%                                                         Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                                         Machanical Engineering, 
%                                                         Energy Science Center,
%                                                         Zürich, July 2016

%--------------------------------------------------------------------------------------------------%
%                                    Fitting of size parameters                                    %
%--------------------------------------------------------------------------------------------------%

function res = MGTsize_fittingPWA (Q, P, S, a, idx, k)

% INITIALIZATION
% k(6) = 1;
% k(7) = 0;
res  = 0;

% SOLVING
for i = 1 : length(S)
  j = i + idx;
  alpha  = fieldnames(a);
  alpha1 = char(alpha(j));
  alpha2 = extractfield(a, alpha1)';
  x      = floor(length(alpha2)/2);
  y      = ceil(length(alpha2)/2);
  alpha3 = alpha2(1:x);
  alpha4 = alpha2(y:end);
  beta   = fieldnames(Q);
  beta1  = char(beta(j));
  beta2  = extractfield(Q, beta1)';  
  gamma  = fieldnames(P);
  gamma1 = char(gamma(j));
  gamma2 = extractfield(P, gamma1)'; 
  q1     = (k(4) * S(i) + k(5)) * (1 - alpha3) + (k(6) * S(i) + k(7)) * alpha3;
  p1     = k(1) * q1 + k(2) * S(i) + k(3);
  q2     = (k(4) * S(i) + k(5)) * (1 - alpha4) + (k(6) * S(i) + k(7)) * alpha4;
  p2     = k(8) * q2 + k(9) * S(i) + k(10);
  q      = [q1; q2];
  p      = [p1; p2];

  res1 = 1/length(beta2) * ((beta2 - q).^2 / S(i));
  res2 = 1/length(gamma2) * ((gamma2 - p).^2 / S(i));
  res  = res + sum(res1 + res2);
  
end

end

%%                                    CODING HYSTORY

% ((k(4) * S(i) + k(5)) * (1 - alpha3) + (k(6) * S(i) + k(7)) * alpha3)