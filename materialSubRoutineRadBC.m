function [eps, deps_u] = materialSubRoutineRadBC(u, params)
%-----------------------------------------------------------------------------------
% Description: This function calculates the radiative boundary material properties and derivatives.
%               
% Input Variables : u = temeprature on integration point.
%                   params = struct containing material parameters.
%
%
% Output Variables : eps = radiative coefficient.
%                    deps_u = derivative of radiative coefficient wrt to temperature.   
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Emissivety coeff parameters
% eps_1 = params.eps_1;
% eps_2 = params.eps_2;

%% Calculate integration point quantities
% Emissivety coeff and concentration derivative
% eps = eps_2+eps_1.*u;
% deps_u = eps_1;


%% Calculate table quantities
GFRP2 = params.GFRP2.GFRP2;
[eps, deps_u] = getTableValues(GFRP2.epsilon, u);

end

