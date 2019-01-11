function [igMatParams] = materialSubRoutineConvBC(u, u_a)
%-----------------------------------------------------------------------------------
% Description: This function calculates the convective boundary material properties and derivatives. 
%               
% Input Variables : u = temeprature on integration point.
%                   u_a = ambient temperature on boundary.
%
% Output Variables : matParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Convective coeff parameters
k_g = 0.03;
Pr = 0.722;
g = 9.81;
beta = 3.43e-3;
nu = 1.57e-5;

%% Calculate integration point quantities

% Convective coeff and concentration derivative
h = 0.14*k_g*(Pr*(g*beta)/nu*(u - u_a))^(1/3);
dh_u = (0.46667*g*beta*Pr*k_g)/((nu*(g*beta*Pr*(u-u_a))/nu)^(2/3));

% Assign calculated quantities to struct
igMatParams.h = h;
igMatParams.dh_u = dh_u;

end

