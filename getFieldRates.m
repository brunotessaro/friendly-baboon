function [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u, c, params)
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
A = params.A;
eta = params.eta;

% Initialize global matrices
Fu = zeros(nNds,1);
Fc = zeros(nNds,1);
Ku = zeros(nNds);
Kc = zeros(nNds);
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
    Kc_e = zeros(size(Te,2));
    Fu_e = zeros(size(Te,2),1);
    Fc_e = zeros(size(Te,2),1);
    
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
            [igMatParams] = materialSubRoutine(psi.sf(ig,:)*u(Te), psi.sf(ig,:)*c(Te));
            k = igMatParams.k;
            Cp = igMatParams.Cp;
            rho = igMatParams.rho;
            omega = igMatParams.omega;      
            
            % Calculate temperature eq. elemental matrices
            Ku_e = Ku_e + gWts(ig)*(M*rho*Cp)*det(J);
            Fu_e = Fu_e + gWts(ig)*(-B'*k*B*u(Te) + M*Q(Te))*det(J);            

            % Calculate concentration eq. elemental matrices
            Kc_e = Kc_e + gWts(ig)*M*det(J);
            Fc_e = Fc_e + gWts(ig)*(psi.sf(ig,:)'*A*(1-psi.sf(ig,:)*c(Te))^eta*omega)*det(J); 
    
    end 
    
    % Assembly of global matrices
    Fu(Te) = Fu(Te) + Fu_e;
    Fc(Te) = Fc(Te) + Fc_e;
    Ku(Te,Te) = Ku(Te, Te) + Ku_e;
    Kc(Te,Te) = Kc(Te, Te) + Kc_e;
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
                q_in = bcInfo{i}{3};
                
                % Calculate temperature eq. elemental matrices for flux boundary
                Fu_e = Fu_e + gWts(ig)*(psi.sf(ig,:)'*q_in)*det(J);
            
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
                
                % Get ambient temperature and convection coeff
                [igMatParams] = materialSubRoutine(psi.sf(ig,:)*u(Te), psi.sf(ig,:)*c(Te));
                h = igMatParams.h;
                u_a = bcInfo{i}{3};
                
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
                
                % Get ambient temperature and emissivity coeff
                [igMatParams] = materialSubRoutine(psi.sf(ig,:)*u(Te), psi.sf(ig,:)*c(Te));
                eps = igMatParams.eps;
                u_a = bcInfo{i}{3};           
                
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
udot(nodeInfo.free) = Ku(nodeInfo.free,nodeInfo.free)\Fu(nodeInfo.free);
cdot = Kc\Fc;

% Assign to struct
rates.udot = udot;
rates.cdot = cdot;

end



