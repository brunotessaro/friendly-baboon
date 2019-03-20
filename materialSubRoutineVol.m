function [igMatParams] = materialSubRoutineVol(u, params)
%-----------------------------------------------------------------------------------
% Description: This function calculates the volumetric material properties and derivatives. 
%               
% Input Variables : u = temeprature on integration point.
%                   params = struct containing material parameters.
%
% Output Variables : igMatParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 
GFRP = params.GFRP.GFRP;

%% Calculate integration point quantities
[k, dk_u] = getTableValues(GFRP.kx, u);
[Cp, dCp_u] = getTableValues(GFRP.cp, u);
[rho, drho_u] = getTableValues(GFRP.rho, u);

% Assign calculated quantities to struct
igMatParams.k = k;
igMatParams.Cp = Cp;
igMatParams.rho = rho;
igMatParams.dk_u = dk_u;
igMatParams.dCp_u = dCp_u;
igMatParams.drho_u = drho_u;


end

