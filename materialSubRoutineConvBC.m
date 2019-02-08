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
h = 0.14E0.*k_g.*(beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a)).^(1/3);
h = 0.14E0.*k_g.*abs((beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a)).^(1/3))*sign((beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a)));

if (u_a == u)
    dh_u = 0;
else
    dh_u = 0.466667E-1.*beta.*g.*k_g.*nu.^(-1).*Pr.*(beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a)).^(-2/3);
    dh_u = 0.466667E-1.*beta.*g.*k_g.*nu.^(-1).*Pr.*(abs((beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a)).^(1/3))*sign((beta.*g.*nu.^(-1).*Pr.*(u+(-1).*u_a))))^(-2);

    
end

end

