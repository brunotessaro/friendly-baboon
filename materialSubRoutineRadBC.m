function [igMatParams] = materialSubRoutineRadBC(u, params)
%-----------------------------------------------------------------------------------
% Description: This function calculates the radiative boundary material properties and derivatives.
%               
% Input Variables : u = temeprature on integration point.
%                   params = struct containing material parameters.
%
%
% Output Variables : matParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Emissivety coeff parameters
eps_1 = params.eps_1;
eps_2 = params.eps_2;

%% Calculate integration point quantities

% Emossivety coeff and concentration derivative
eps = eps_1*u+eps_2;
deps_u = eps_1;

% Assign calculated quantities to struct
igMatParams.eps = eps;
igMatParams.deps_u = deps_u;

end

