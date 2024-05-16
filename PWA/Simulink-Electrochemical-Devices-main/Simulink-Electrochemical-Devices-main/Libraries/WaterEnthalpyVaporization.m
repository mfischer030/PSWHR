%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                     Water enthalpy of vaporization                      %
%-------------------------------------------------------------------------%

function hH2Ovap = WaterEnthalpyVaporization (p_WHSG, T)

if p_WHSG > 0

 % Water enthalpy of vaporization (p) [J / mol]
 psi1     = -101.035293;        % Fitting coefficient - Thermodynamic
 psi2     = 707.851918;         % Fitting coefficient - Thermodynamic
 psi3     = -2426.187284;       % Fitting coefficient - Thermodynamic
 psi4     = +42478.862098;      % Fitting coefficient - Thermodynamic
 hH2Ovap = (psi1 * p_WHSG^3 + psi2 * p_WHSG^2 + psi3 * p_WHSG + psi4);

else
    
 % Water enthalpy of vaporization (T) [J / mol]
 psi1     = - 5.115e-13;        % Fitting coefficient - Thermodynamic
 psi2     = 1.34605e-09;        % Fitting coefficient - Thermodynamic
 psi3     = - 1.45420542e-06;   % Fitting coefficient - Thermodynamic
 psi4     = 8.2503854758e-04;   % Fitting coefficient - Thermodynamic
 psi5     = - 0.259106785469;   % Fitting coefficient - Thermodynamic
 psi6     = 42.6523745786349;   % Fitting coefficient - Thermodynamic
 psi7     = - 2826.638120809;   % Fitting coefficient - Thermodynamic
 hH2Ovap = (psi1 * T^6 + psi2 * T^5 + psi3 * T^4 + psi4 * T^3 + ...
    psi5 * T^2 + psi6 * T + psi7);

end

end