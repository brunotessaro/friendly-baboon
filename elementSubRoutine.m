function [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u, u_n, rates, u_inf, params)
%-----------------------------------------------------------------------------------
% Description: This function creates the residual and tangent matrix of the problem. 
%               
% Input Variables : nodeInfo = struct containing node related info.
%                   elemInfo = struct containing element related info.
%                   bcInfo = struct containing boundary conditions info.
%                   u, u_n = temperature in current and previous time step.
%                   rates = struct containing the field rates at the previous time step.
%                   u_inf = input fire temperature.
%                   params = struct containing information of physical and numerical parameters.
%
% Output Variables : r = residual vector.
%                    K = tangent matrix.
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Get mesh parameters
X = nodeInfo.X;
T = elemInfo.T;
elemType = elemInfo.elemType;
volElemIdx = elemInfo.volElemIdx;
volIntRule = elemInfo.volIntRule;
nNds = size(X,1);

% Get physical and numerical parameters
sig = params.sig;
dt = params.dt;
alpha = params.alpha;

% Get rates at previous time step
udot_n = rates.udot;

% Initialize global tangents and resdiuals
r_u = zeros(nNds,1);
K_uu = zeros(nNds);

%% Assembly of volume elements contributions

% Loop trough volume elements
for i=1:size(volElemIdx,1)
    
    % Retrieve element data
    iElem = volElemIdx(i);
    [psi, gWts, nElNds, dim] = elementCall(elemType(iElem), volIntRule);
    Te = T(iElem,1:nElNds);
    Xe = X(Te,:);
    
    % Initialize elemental matrices
    re_u = zeros(size(Te,2),1);
    Ke_uu = zeros(size(Te,2));
    
    % Loop trough integration points
    for ig=1:size(gWts,2)
        
        % Get jacobian matrix
        J = jacobCalc(elemType(iElem), ig, Xe, psi);
        
        % Get auxiliary matrices
        if dim == 1
            B = J\psi.dcsi(ig,:);
        elseif dim == 2
            B = J\[psi.dcsi(ig,:) ; psi.deta(ig,:)];
        end
        M = psi.sf(ig,:)'*psi.sf(ig,:);
        
        % Get material point quantities
        [igMatParams] = materialSubRoutineVol(psi.sf(ig,:)*u(Te), params);
        k = igMatParams.k;
        Cp = igMatParams.Cp;
        rho = igMatParams.rho;
        dk_u = igMatParams.dk_u;
        dCp_u = igMatParams.dCp_u;
        drho_u = igMatParams.drho_u;
        
        % Calculate residual vectors for temperature and concentration
        re_u = re_u + gWts(ig)*(M*rho*Cp*(u(Te) - u_n(Te) - dt*(1-alpha)*udot_n(Te)) ...
                                + alpha*dt*(B'*k*B*u(Te)))*det(J);
        
        % Calculate tangent matrices
        Ke_uu = Ke_uu + gWts(ig)*(M*rho*Cp ...
                                  + M*(drho_u*Cp + dCp_u*rho)*(psi.sf(ig,:)*(u(Te) - u_n(Te) - dt*(1-alpha)*udot_n(Te))) ...
                                  + alpha*dt*(B'*k*B + (B'*dk_u*(B*u(Te))*psi.sf(ig,:))))*det(J);
                                      
    end
    
    % Assembly of global matrices
    r_u(Te) = r_u(Te) + re_u;
    K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;

end

%% Assembly of boundary element contributions

% Loop trough boundary physical tags
for i=1:size(bcInfo,2)
    
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    if bcInfo{i}{2} == 1
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{6},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{6}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), bcInfo{i}{5});
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices
            re_u = zeros(size(Te,2),1);
            
            % Loop trough integration points
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix 
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get normal flux
                qb = bcInfo{i}{3};
                
                % Calculate temperature eq. residual for flux boundary
                re_u = re_u + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*qb)*det(J);
                
            end
            
            % Assembly of global matrices
            r_u(Te) = r_u(Te) + re_u;
        end
        
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    elseif bcInfo{i}{2} == 2
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{6}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), bcInfo{i}{5});
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices         
            re_u = zeros(size(Te,2),1);
            Ke_uu = zeros(size(Te,2));
            
            % Loop trough integration points
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get ambient temperature and convection coeff, this is problem dependent and the if's
                % have to be hard coded for different kinds of boundary conditions.
                if bcInfo{i}{1} == 10
                    u_a = u_inf;
                    h = bcInfo{i}{4};
                    dh_u = 0;
                elseif bcInfo{i}{1} == 12
                    u_a = bcInfo{i}{3};
                    [h, dh_u] = materialSubRoutineConvBC(psi.sf(ig,:)*u(Te), u_a, params);
                else
                    u_a = bcInfo{i}{3};
                    h = bcInfo{i}{4};
                    dh_u = 0;
                end
                
                % Calculate temperature eq. residual for convective bc
                re_u = re_u + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*h*(u_a - psi.sf(ig,:)*u(Te)))*det(J);
                
                % Calculate tangents
                Ke_uu = Ke_uu + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*(dh_u*(u_a - psi.sf(ig,:)*u(Te)) - h)*psi.sf(ig,:))*det(J);
                
            end
            
            % Assembly of global matrices
            r_u(Te) = r_u(Te) + re_u;
            K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;
        end
        
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    elseif bcInfo{i}{2} == 3
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{6}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), bcInfo{i}{5});
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices         
            re_u = zeros(size(Te,2),1);
            Ke_uu = zeros(size(Te,2));
            
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix  
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get ambient temperature and emissivity coeff, this is problem dependent and the if's
                % have to be hard coded for different kinds of boundary conditions.
                if bcInfo{i}{1} == 11
                    [eps, deps_u] = materialSubRoutineRadBC(psi.sf(ig,:)*u(Te), params);
                    u_a = u_inf;
                elseif bcInfo{i}{1} == 13
                    [eps, deps_u] = materialSubRoutineRadBC(psi.sf(ig,:)*u(Te), params);
                    u_a = bcInfo{i}{3};
                else
                    u_a = bcInfo{i}{3};
                    eps = bcInfo{i}{4};
                    deps_u = 0;
                end
                
                % Calculate temperature eq. residual for radiative bc and tangets
                re_u = re_u + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*sig*eps*((u_a + 273.15)^4 - (psi.sf(ig,:)*(u(Te) + 273.15))^4))*det(J);
                
                % Calculate tangents
                Ke_uu = Ke_uu + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*sig*(deps_u*((u_a+273.15)^4 - (psi.sf(ig,:)*(u(Te) + 273.15))^4)  ...
                                          -4*eps*(psi.sf(ig,:)*(u(Te) + 273.15))^3)*psi.sf(ig,:))*det(J);

            end
            
            % Assembly of global matrices
            r_u(Te) = r_u(Te) + re_u;
            K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;
        end
    end
end

%% Grouping global matrices for newton iteration
r = r_u;
K = K_uu;

end

