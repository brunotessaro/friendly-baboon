function [igMatParams] = materialSubRoutine(u, c)
%-----------------------------------------------------------------------------------
% Description: This function creates the residual and tangent matrix of the problem. 
%               
% Input Variables : u = temeprature on integration point.
%                   c = concentration on integration point.
%
% Output Variables : matParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% For later

%% Calculate integration point quantities
k = 1;
Cp = 1;
rho = 1;
h = 1;
eps = 1;
omega = 1;

% Assign calculated quantities to struct
igMatParams.k = k;
igMatParams.Cp = Cp;
igMatParams.rho = rho;
igMatParams.h = h;
igMatParams.eps = eps;
igMatParams.omega = omega;


end

