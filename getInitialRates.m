function [dudt] = getInitialRates(nodeInfo, elemInfo, bcInfo, u, params)
%% Get auxiliary paramenters 

% Get mesh parameters
X = nodeInfo.X;
T = elemInfo.T;
elemType = elemInfo.elemType;
volElemIdx = elemInfo.volElemIdx;

% Initialize sizes
nNds = size(X,1);

% Get physical parameters
k = params.k;
Q = params.Q;
sig = params.sig;
rho = params.rho;
Cp = params.Cp;

% Initialize global matrices
A = zeros(nNds,1);
M = zeros(nNds);
dudt = zeros(nNds,1);

%% Assembly of volume elements contributions

% Loop trough volume elements
for i=1:size(volElemIdx,1)
    
    % Retrieve element data and initialize elemental matrices
    iElem = volElemIdx(i);
    [psi, gWts, nElNds, dim] = elementCall(elemType(iElem));
    Te = T(iElem,1:nElNds);
    Xe = X(Te,:);
    Me = zeros(size(Te,2));
    Ae = zeros(size(Te,2),1);
    
    % Loop trough element gauss points
    for ig=1:size(gWts,2)
            J = jacobCalc(elemType(iElem), ig, Xe, psi);
            
            % Assemble B matrix depending on dimension.
            if dim == 1
                B = J\psi.xi(ig,:); 
            elseif dim == 2
                B = J\[psi.xi(ig,:) ; psi.eta(ig,:)];
            end
            
            % Calculate elemental mass matrix and force vector
            Me = Me + gWts(ig)*(psi.bf(ig,:)'*rho*Cp*psi.bf(ig,:))*det(J);
            Ae = Ae + gWts(ig)*(-B'*k*B*u(Te) + psi.bf(ig,:)'*(psi.bf(ig,:)*Q(Te)))*det(J);            

    end 
    
    % Assembly of mass matrix and force vector
    A(Te) = A(Te) + Ae;
    M(Te,Te) = M(Te, Te) + Me;
end

%% Assembly of boundary elements contributions

% Loop trough boundary physical tags
for i=1:size(bcInfo,2)
    
    % Check the kind of boundary condition
    if bcInfo{i}{2} == 1
        
        % Loop trough elements in boundary tag
        for j=1:size(bcInfo{i}{1},1)
            
            % Retrieve element data and initialize elemental matrices
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem));
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            Ae = zeros(size(Te,2),1);
            
            % Loop trough element gauss points
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                Ae = Ae + gWts(ig)*(psi.bf(ig,:)'*bcInfo{i}{3})*det(J);
            end
        	% Assembly of residual and tangent matrix
            A(Te) = A(Te) + Ae;
        end
        
    % Same procedure for convective boundary  
    elseif bcInfo{i}{2} == 2
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem));
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            Ae = zeros(size(Te,2),1);
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                Ae = Ae + gWts(ig)*(psi.bf(ig,:)'*bcInfo{i}{4}*(bcInfo{i}{3} - psi.bf(ig,:)*u(Te)))*det(J);
            end 
            A(Te) = A(Te) + Ae;
        end
        
    % Same procedure for radiative boundary      
    elseif bcInfo{i}{2} == 3
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem));
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            Ae = zeros(size(Te,2),1);
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                Ae = Ae + gWts(ig)*(psi.bf(ig,:)'*sig*bcInfo{i}{4}*(bcInfo{i}{3}^4 - (psi.bf(ig,:)*u(Te))^4))*det(J);
            end
            A(Te) = A(Te) + Ae;
        end
    end
end

dudt(nodeInfo.free) = M(nodeInfo.free,nodeInfo.free)\A(nodeInfo.free);


end



