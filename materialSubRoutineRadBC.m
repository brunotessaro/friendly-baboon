function [igMatParams] = materialSubRoutineRadBC(u)
%-----------------------------------------------------------------------------------
% Description: This function calculates the radiative boundary material properties and derivatives.
%               
% Input Variables : u = temeprature on integration point.
%
% Output Variables : matParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Emissivety coeff parameters
eps_a = 2.0408163265306E-4;
eps_b = 0.74591836734694;

%% Calculate integration point quantities

% Emossivety coeff and concentration derivative
eps = eps_a*u+eps_b;
deps_u = eps_a;

% Assign calculated quantities to struct
igMatParams.eps = eps;
igMatParams.deps_u = deps_u;

end

