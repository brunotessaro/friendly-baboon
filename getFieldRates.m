function [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u, u_inf, params)
%-----------------------------------------------------------------------------------
% Description: This function calculate the temperature and concentration rates in a given time step. 
%               
% Input Variables : nodeInfo = struct containing node related info.
%                   elemInfo = struct containing element related info.
%                   bcInfo = struct containing boundary conditions info.
%                   u = temperature on the desired time step.
%                   u_inf = input fire temperature.
%                   params = struct containing information of physical and numerical parameters.
%
% Output Variables : rates = struct containing the rates in the desired time step.
%            
%-----------------------------------------------------------------------------------
%% Get auxiliary paramenters 

% Get mesh parameters
X = nodeInfo.X;
T = elemInfo.T;
elemType = elemInfo.elemType;
volElemIdx = elemInfo.volElemIdx;
nNds = size(X,1);

% Get physical parameters
Q = params.Q;
sig = params.sig;

% Initialize global matrices
Fu = zeros(nNds,1);
Ku = zeros(nNds);
udot = zeros(nNds,1);

%% Assembly of volume element contributions

% Loop trough volume elements
for i=1:size(volElemIdx,1)
    
    % Retrieve element data
    iElem = volElemIdx(i);
    [psi, gWts, nElNds, dim] = elementCall(elemType(iElem), [2 2]);
    Te = T(iElem,1:nElNds);
    Xe = X(Te,:);
    
    % Initialize elemental matrices
    Ku_e = zeros(size(Te,2));
    Fu_e = zeros(size(Te,2),1);
    
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
            
            % Calculate temperature eq. elemental matrices
            Ku_e = Ku_e + gWts(ig)*(M*rho*Cp)*det(J);
            Fu_e = Fu_e + gWts(ig)*(-B'*k*B*u(Te))*det(J);            
    
    end 
    
    % Assembly of global matrices
    Fu(Te) = Fu(Te) + Fu_e;
    Ku(Te,Te) = Ku(Te, Te) + Ku_e;
end

%% Assembly of boundary element contributions

% Loop trough boundary physical tags
for i=1:size(bcInfo,2)
    
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    if bcInfo{i}{2} == 1
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 0]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices
            Fu_e = zeros(size(Te,2),1);
            
            % Loop trough integration points
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix 
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get normal flux
                qb = bcInfo{i}{3};
                
                % Calculate temperature eq. elemental matrices for flux boundary
                Fu_e = Fu_e + gWts(ig)*(psi.sf(ig,:)'*qb)*det(J);
            
            end
            
        	% Assembly of global matrices
            Fu(Te) = Fu(Te) + Fu_e;  
        end
        
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    elseif bcInfo{i}{2} == 2
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 1]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices         
            Fu_e = zeros(size(Te,2),1);
            
             % Loop trough integration points
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix 
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get ambient temperature and convection coeff, this is problem dependent and the if's
                % have to be hard coded for different kinds of boundary conditions.
                if bcInfo{i}{5} == 10
                    u_a = u_inf;
                    h = bcInfo{i}{4};
                elseif bcInfo{i}{5} == 12
                    u_a = bcInfo{i}{3};
                    [h, ~] = materialSubRoutineConvBC(psi.sf(ig,:)*u(Te), u_a, params);
                else
                    u_a = bcInfo{i}{3};
                    h = bcInfo{i}{4};
                end
                
                % Calculate temperature eq. elemental matrices for convective boundary
                Fu_e = Fu_e + gWts(ig)*(psi.sf(ig,:)'*h*(u_a - psi.sf(ig,:)*u(Te)))*det(J);
            end
            
            % Assembly of global matrices
            Fu(Te) = Fu(Te) + Fu_e;
        end
        
    % Check the kind of boundary condition (1: flux, 2: convective, 3: radiative)
    elseif bcInfo{i}{2} == 3
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 4]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            % Initialize elemental matrices         
            Fu_e = zeros(size(Te,2),1);
            
            for ig=1:size(gWts,2)
                
                % Get jacobian matrix                 
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                % Get ambient temperature and emissivity coeff, this is problem dependent and the if's
                % have to be hard coded for different kinds of boundary conditions.
                if bcInfo{i}{5} == 11
                    [eps, ~] = materialSubRoutineRadBC(psi.sf(ig,:)*u(Te), params);
                    u_a = u_inf;
                elseif bcInfo{i}{5} == 13
                    [eps, ~] = materialSubRoutineRadBC(psi.sf(ig,:)*u(Te), params);
                    u_a = bcInfo{i}{3};
                else
                    u_a = bcInfo{i}{3};
                    eps = bcInfo{i}{4};
                end         
                
                % Calculate temperature eq. elemental matrices for radiative boundary
                Fu_e = Fu_e + gWts(ig)*(psi.sf(ig,:)'*sig*eps*(u_a^4 - (psi.sf(ig,:)*u(Te))^4))*det(J);
            end
            
            % Assembly of global matrices
            Fu(Te) = Fu(Te) + Fu_e;
        end
    end
end

%% Calculate rates

% Solve the system to get rates
udot(nodeInfo.free(1:nNds)) = Ku(nodeInfo.free(1:nNds),nodeInfo.free(1:nNds))\Fu(nodeInfo.free(1:nNds));

% Assign to struct
rates.udot = udot;

end



