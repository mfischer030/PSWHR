%-------------------------------------------------------------------------%
%                                 IMES                                    %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                                Process Engineering Institute, 
%                                Energy Science center,
%                                Zürich, March 2016

%-------------------------------------------------------------------------%
%                                 PEMEC                                   %
%-------------------------------------------------------------------------%

%%                       DESCRIPTION AND ASSUMPTIONS

% Thermodynamic model of a proton exchange membrane electrolyzer (PEMEC).
% A lumped, first-principle approach is followed. Static and dynamic 
% behavior are described, and both electrochemical and thermal features of 
% the PEM electrolyzer are captured.

% The data represent a big scale PEM electrolyzer, up to 100 kW.
% It is modeled based on the commercial product Silyzer100 of Siemens
% (100 kW).
% The operative power can range in 5-100% of the rated power.
% The hydrogen and oxygen are assumed to be produced at 40 bar.
% A reference temperature of 70 °C is chosen as the commercial product is
% assumed to be controlled to work close to the design temperature in the
% whole range of power.
% For the moment, heat integration is not considered for the electrolyzer,
% although the device can produce heat while producing hydrogen at high
% power.
% A first order dynamics is used to simulate the time response of both the
% electrical and thermal produced powers. The time constant for the
% electrical power is typically 1 order of magnitude smaller than that of
% the thermal power.
 
% Modeling parameter were chosen from the literature. The following
% assumptions were made:
% 1. 0-D model: greybox approach;
% 2. Ideal gases;
% 3. Negligible pressure drop along the channels;
% 4. Constant operating temperature (design conditions).

%% -----------------------------------------------------------------------%
%                                REFERENCES                               %
%-------------------------------------------------------------------------%

% [1] Oystein Ullberg, Modeling of advanced alkaline electrolysers: a
% system simulation approach, Hydrogen Energy, 2003

% [2] Garcia-Valverde R., Espinosa N., Urbina A., Simple PEM water 
% electrolyser model and experimental validation, International Journal 
% of Hydrogen Energy, Volume 37, Issue 2, 2012.

% [3] Marangio F., Santarelli M., Cali M., Simple Theoretical model and 
% experimental analysis of a high pressure PEM water electrolyser for 
% hydrogen production, International Journal of Hydrogen Energy, 2009.

% [4] H. Kim, M. Park, and L.S. Lee, One dimensional dynamic modeling of a
% high-pressure water electrolysis system for hydrogen production, Int. J.
% Hydrogen energy, vol. 38, pp. 2596-2609, 2013.
