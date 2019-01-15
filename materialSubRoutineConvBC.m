function [h, dh_u] = materialSubRoutineConvBC(u, u_a, params)
%-----------------------------------------------------------------------------------
% Description: This function calculates the convective boundary material properties and derivatives. 
%               
% Input Variables : u = temeprature on integration point.
%                   u_a = ambient temperature on boundary.
%                   params = struct containing material parameters.
%
% Output Variables : h = convection coefficient.
%                    dh_u = derivative of convection coefficient wrt to temperature.
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Convective coeff parameters
k_g = params.k_g;
Pr = params.Pr;
g = params.g;
beta = params.beta;
nu = params.nu;

%% Calculate integration point quantities

% Convective coeff and concentration derivative
h = 0.14*k_g*(Pr*(g*beta)/nu*(u - u_a))^(1/3);
dh_u = (0.46667*g*beta*Pr*k_g)/((nu*(g*beta*Pr*(u-u_a))/nu)^(2/3));

end

