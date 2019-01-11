function [igMatParams] = materialSubRoutineVol(u, c)
%-----------------------------------------------------------------------------------
% Description: This function calculates the volumetric material properties and derivatives. 
%               
% Input Variables : u = temeprature on integration point.
%                   c = concentration on integration point.
%
% Output Variables : matParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Conduction coeff parameters
k_a = 0.1;
k_b = 0.35;

% Rho parameters
rho_a = 1141;
rho_b = 1870;

% Specific heat capacity parameters
Cp_a = 840;
Cp_b = 1170;
Cp_d = 234e3;

% Omega parameters
R = 8.314;
E = 77878;

%% Calculate integration point quantities

% Conduction coeff and concentration derivative
k = (k_a*k_b)/(k_a*(1-c)+k_b*c)*eye(2);
dk_c = (k_a*k_b*(k_b - k_a))/((1-c)*(k_a + c*k_b))^2;

% Specific heat and concentration derivative
f = (Mi*(1 - c))/(Mi*(1 - c) + Mf*c);

Cp = Cp_b*f + Cp_a(1-f) + (c - c_n)/(u - u_n)*Cp_d;

dCp_c = ((Cp_a - Cp_b)*Mi*Mf)/(c*(Mf - Mi) + Mi)^2+Cp_d/(u - u_n);
dCp_u = ((c_n-c)*Cp_d)/(u-u_n^2);

% Density and concentration derivative
rho = (1-c)*rho_b + c*rho_a;
drho_c = rho_a - rho_b;

% Omega and temperature derivative
omega = exp(-E/(R*u));
domega_u = (E*exp(-E/(R*u)))/(R*u^2);

% Assign calculated quantities to struct
igMatParams.k = k;
igMatParams.Cp = Cp;
igMatParams.rho = rho;
igMatParams.omega = omega;
igMatParams.dk_c = dk_c;
igMatParams.dCp_c = dCp_c;
igMatParams.dCp_u = dCp_u;
igMatParams.drho_c = drho_c;
igMatParams.domega_u = domega_u;

end

