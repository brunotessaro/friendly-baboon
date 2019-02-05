function [igMatParams] = materialSubRoutineVol(u, u_n, c, c_n, params)
%-----------------------------------------------------------------------------------
% Description: This function calculates the volumetric material properties and derivatives. 
%               
% Input Variables : u = temeprature on integration point.
%                   u_n = temeprature on integration point on previous time step.
%                   c = concentration on integration point.
%                   c_n = concentration on integration point on previous time step.
%                   params = struct containing material parameters.
%
% Output Variables : igMatParams = struct containing all the quantities calculated in the gauss pt.
%                    
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

k_a = params.k_a;
k_b = params.k_b;
rho_a = params.rho_a;
rho_b = params.rho_b;
Cp_a = params.Cp_a;
Cp_b = params.Cp_b;
Cp_d = params.Cp_d;

%% Calculate integration point quantities

% The quantities here presented are calculated symbolic on Mathematica and exported aumtomaticly to
% avoid typing or derivative errors.

% Conduction coeff and concentration derivative
k = (c.*k_a.^(-1)+(1+(-1).*c).*k_b.^(-1)).^(-1);
dk_c = k_a.*(k_a+(-1).*k_b).*k_b.*(k_a+(-1).*c.*k_a+c.*k_b).^(-2);

% Density and concentration derivative
rho = c.*rho_a+(1+(-1).*c).*rho_b;
drho_c = rho_a+(-1).*rho_b;

% Specific heat and concentration derivative
% Obs: since the author uses a numerical trick to calculate the derivative dc/du one has to place
% countermeasures in the case u does not change (u == u_n) in order to avoid division by 0.
if (u == u_n) 
    Cp = (1+(-1).*c).*Cp_b.*rho_b.*(c.*rho_a+(1+(-1).*c).*rho_b).^(-1)+Cp_a.*(1+(-1).* ...
         (1+(-1).*c).*rho_b.*(c.*rho_a+(1+(-1).*c).*rho_b).^(-1));
    dCp_c = (Cp_a+(-1).*Cp_b).*rho_a.*rho_b.*(c.*(rho_a+(-1).*rho_b)+rho_b).^(-2);
    dCp_u = 0;
else
    Cp = (1+(-1).*c).*Cp_b.*rho_b.*(c.*rho_a+(1+(-1).*c).*rho_b).^(-1)+Cp_a.*(1+(-1).* ...
         (1+(-1).*c).*rho_b.*(c.*rho_a+(1+(-1).*c).*rho_b).^(-1))+(c+(-1).*c_n).* ...
         Cp_d.*(u+(-1).*u_n).^(-1);   
    dCp_c = (Cp_a+(-1).*Cp_b).*rho_a.*rho_b.*(c.*(rho_a+(-1).*rho_b)+rho_b).^(-2)+Cp_d.*(u+(-1).*u_n).^(-1);
    dCp_u = ((-1).*c+c_n).*Cp_d.*(u+(-1).*u_n).^(-2);
end

% Phase-field parameters

% Symetric part of potential g(c)
% g = c^2*(1-c)^2;
dg_c = 2*c*(2*c^2 - 3*c + 1);
dg_cc = 3*(2*c-1)^2 -1;

% Tilted part of potential g(c)

% Type I z(c)
% z = c^3*(10 - 15*c + 6*c^2);
% dz_c = 30*c^2*(c-1)^2;
% dz_cc = 60*c*(2*c^2 - 3*c + 1);

% Type II z(c)
% z = c^2(3-2c);
dz_c = -6*c*(c-1);
dz_cc = 6-12*c;

% f function
% ------------------------------
% cte = -0.0000001;
% f = cte*u;
% df_u = cte;
cte = -0.163934426229508;
f = cte;
df_u = 0;

% Assign calculated quantities to struct
igMatParams.k = k;
igMatParams.dk_c = dk_c;
igMatParams.Cp = Cp;
igMatParams.dCp_c = dCp_c;
igMatParams.dCp_u = dCp_u;
igMatParams.rho = rho;
igMatParams.drho_c = drho_c;

igMatParams.dg_c = dg_c;
igMatParams.dg_cc = dg_cc;
igMatParams.dz_c = dz_c;
igMatParams.dz_cc = dz_cc;
igMatParams.f = f;
igMatParams.df_u = df_u;

end

