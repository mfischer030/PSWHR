%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                        Water saturation pressure                        %
%-------------------------------------------------------------------------%

function pH2Osat = WaterSaturationPressure (T)

chi1     = 0.00129697;   % Fitting coefficient - Thermodynamic
chi2     = - 1.529053;   % Fitting coefficient - Thermodynamic
chi3     = 681.731481;   % Fitting coefficient - Thermodynamic
chi4     = - 136025.0;   % Fitting coefficient - Thermodynamic
chi5     = 10234070.7;   % Fitting coefficient - Thermodynamic

% Water saturation pressure [atm]
pH2Osat = (chi1 * T^4 + chi2 * T^3 + chi3 * T^2 + chi4 * T + chi5) / ...
    101325;

end