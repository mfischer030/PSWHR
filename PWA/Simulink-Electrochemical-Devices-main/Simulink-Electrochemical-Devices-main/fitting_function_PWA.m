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

function res = fitting_function(IN, OUT, S, a, N_bp, k)

% INITIALIZATION
res = 0;

% LINEAR FITTING
if N_bp < 2

  for i = 1 : length(S) 
    in  = (k(4) * S(i) + k(5)) .* (1 - a{i}) + (k(6) * S(i) + k(7)) .* a{i};
    out = k(1) * in + k(2) * S(i) + k(3);

    res1 = 1/length(IN{i}) * ((IN{i} - in).^2 / S(i));
    res2 = 1/length(OUT{i}) * ((OUT{i} - out).^2 / S(i));
    res  = res + sum(res1 + res2);
  end
  
% PWA FITTING
else
  theta = cell(1, N_bp - 1);
  in    = cell(1, N_bp - 1);
  out   = cell(1, N_bp - 1);
  for i = 1 : length(S) 
    index = round(linspace(1, length(a{i}), N_bp));
    for j = 1 : N_bp - 2
      theta{j} = a{i}(1:index(j+1));
      in{j}    = (k(4+6+j) * S(i) + k(5+6+j)) .* (1 - theta{j}) + (k(6+6+j) * S(i) + k(7+6+j)) .* theta{j};
      out{j}   = k(1+6+j) * in{j} + k(2+6+j) * S(i) + k(3+6+j);
      res1 = 1/length(theta{j}) * ((IN{i}(1:index(j+1)) - in{j}).^2 / S(i));
      res2 = 1/length(theta{j}) * ((OUT{i}(1:index(j+1)) - out{j}).^2 / S(i));
      res  = res + sum(res1 + res2);
      
      if j == N_bp - 2
        theta{j} = a{i}(index(j+1)+1:end);
        in{j}    = (k(4+6+j) * S(i) + k(5+6+j)) .* (1 - theta{j}) + (k(6+6+j) * S(i) + k(7+6+j)) .* theta{j};
        out{j}   = k(1+6+j) * in{j} + k(2+6+j) * S(i) + k(3+6+j);
        res1 = 1/length(theta{j}) * ((IN{i}(index(j+1)+1:end) - in{j}).^2 / S(i));
        res2 = 1/length(theta{j}) * ((OUT{i}(index(j+1)+1:end) - out{j}).^2 / S(i));
        res  = res + sum(res1 + res2);
      end
      
    end
  end
  
end

end