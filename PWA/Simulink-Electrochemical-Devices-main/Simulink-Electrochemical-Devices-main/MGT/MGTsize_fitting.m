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

function res = MGTsize_fitting (Q, P, S, a, idx, k)

% INITIALIZATION
% k(6) = 1;
% k(7) = 0;
k(8)  = 0;
k(9)  = 0;
k(10) = 0;
res   = 0;

% SOLVING
for i = 1 : length(S)
  j = i + idx;
  alpha  = fieldnames(a);
  alpha1 = char(alpha(j));
  alpha2 = extractfield(a, alpha1)';
  beta   = fieldnames(Q);
  beta1  = char(beta(j));
  beta2  = extractfield(Q, beta1)';  
  gamma  = fieldnames(P);
  gamma1 = char(gamma(j));
  gamma2 = extractfield(P, gamma1)'; 
  q      = (k(4) * S(i) + k(5)) * (1 - alpha2) + (k(6) * S(i) + k(7)) * alpha2;
  p      = k(1) * ((k(4) * S(i) + k(5)) .* (1 - alpha2) + (k(6) * S(i) + k(7)) .* alpha2) + k(2) * S(i) + k(3);

  res1 = 1/length(beta2) * ((beta2 - q).^2 / S(i));
  res2 = 1/length(gamma2) * ((gamma2 - p).^2 / (S(i) * 0.28));
  res  = res + sum(res1 + res2);
  
%   for z = 1 : length(p)
%     if p(z) < 0.04
%       res = 1000;
%     end
%   end
  
end

end