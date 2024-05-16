%-------------------------------------------------------------------------%
%                    ALKEC - Paolo Gabrielli - ETH Zurich                 %
%-------------------------------------------------------------------------%

%                                Paolo Gabrielli (gapaolo@ipe.mavt.ethz.ch)
%                              Process Engineering Institute, November 2015


%-------------------------------------------------------------------------%
%        Preliminary ALKEC Design - Hydrogenics HyStat10 - 100 kW         %
%-------------------------------------------------------------------------%

%%                      DESCRIPTION AND ASSUMPTIONS

% Preliminary design of an alkaline electrolytic cell (ALKEC). The cell is 
% modeled as a 0-D grey box component, where physical laws are combined 
% with data from literature. Static and dynamic behavior are described. 
% Electrochemical and thermal features of the alkaline electrolyzer are 
% captured. Modeling parameter were chosen from the literature, by assuming 
% the following:
% 1. 0-D model: greybox approach;
% 2. Ideal gases;
% 3. Negligible pressure drop along the channels;
% 4. Constant operating temperature (design conditions).

%% -----------------------------------------------------------------------%
%                                REFERENCES                               %
%-------------------------------------------------------------------------%

% [1] Oystein Ullberg, Modeling of advanced alkaline electrolysers: a
% system simulation approach, Hydrogen Energy, 2003

% [2] Padulles J., Ault G.W., McDonald, J.R., An integrated SOFC plant
% dynamic model for power system simulations, 2000