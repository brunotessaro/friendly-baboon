function [igMatParams] = materialSubRoutine(u, c, params)
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
k_a = params.k_a;
k_b = params.k_b;

rho_a = params.rho_a;
rho_b = params.rho_b;


params.k_g = 1;


%% Calculate integration point quantities
k = (k_a*k_b)/(k_a*(1-c)+k_b*c);
Cp = 1;
rho = (1-c)*rho_b + c*rho_a;
h = 0.14*k_g*(Pr*(g*beta)/nu*(u - u_a))^(1/3);
eps = 2.0408163265306E-4*u+0.74591836734694;
omega = 1;

dk_c = 1;
dCp_c = 1;
drho_c = 1;
domega_u = 1;
dh_u = 1;
deps_u = 1;

% Assign calculated quantities to struct
igMatParams.k = k;
igMatParams.Cp = Cp;
igMatParams.rho = rho;
igMatParams.h = h;
igMatParams.eps = eps;
igMatParams.omega = omega;

igMatParams.dk_c = dk_c;
igMatParams.dCp_c = dCp_c;
igMatParams.drho_c = drho_c;
igMatParams.domega_u = domega_u;
igMatParams.dh_u = dh_u;
igMatParams.deps_u = deps_u;


end

