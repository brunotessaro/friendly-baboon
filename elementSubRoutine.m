function [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u, u_n, c, c_n, params)
%-----------------------------------------------------------------------------------
% Description: This function creates the residual and tangent matrix of the problem. 
%               
% Input Variables : nodeInfo = struct containing node related info.
%                   elemInfo = struct containing element related info.
%                   bcInfo = struct containing boundary conditions info.
%                   u, u_n = temperature in current and previous time step.
%                   c, c_n = concentration in current and previous time step.
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

% Initialize sizes
nNds = size(X,1);

% Get physical parameters
k = params.k;
Q = params.Q;
sig = params.sig;
rho = params.rho;
Cp = params.Cp;
dt = params.dt;
alpha = params.alpha;

% Initialize global matrices
K = zeros(nNds);
r = zeros(nNds,1);

%% Assembly of volume elements contributions

% Loop trough volume elements
for i=1:size(volElemIdx,1)
    
    % Retrieve element data and initialize elemental matrices
    iElem = volElemIdx(i);
    [psi, gWts, nElNds, dim] = elementCall(elemType(iElem), [2 2]);
    
    Te = T(iElem,1:nElNds);
    Xe = X(Te,:);
    Ke = zeros(size(Te,2));
    re = zeros(size(Te,2),1);
    
    % Loop trough element gauss points
    for ig=1:size(gWts,2)
            J = jacobCalc(elemType(iElem), ig, Xe, psi);
            
            % Assemble B matrix depending on dimension.
            if dim == 1
                B = J\psi.dcsi(ig,:); 
            elseif dim == 2
                B = J\[psi.dcsi(ig,:) ; psi.deta(ig,:)];
            end
            
            M = psi.sf(ig,:)'*rho*Cp*psi.sf(ig,:);
            Kappa = B'*k*B;
            F_Q = M*Q(Te);
            
            re = re + gWts(ig)*(M*(u(Te) - u_n(Te)) ...
                                + dt*(1 - alpha)*(Kappa*u_n(Te) - F_Q) ...
                                + alpha*dt*(Kappa*u(Te) - F_Q))*det(J);     
            
            Ke = Ke + gWts(ig)*(M + alpha*dt*Kappa)*det(J);
            
    end 
    
    % Assembly of residual and tangent matrix
     r(Te) = r(Te) + re;
     K(Te,Te) = K(Te, Te) + Ke;
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
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 0]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            re = zeros(size(Te,2),1);
            
            % Loop trough element gauss points
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                re = re + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*bcInfo{i}{3} ...
                                    -dt*(1-alpha)*psi.sf(ig,:)'*bcInfo{i}{3})*det(J);
            end
        	% Assembly of residual and tangent matrix
            r(Te) = r(Te) + re;
        end
        
    % Same procedure for convective boundary  
    elseif bcInfo{i}{2} == 2
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 1]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            Ke = zeros(size(Te,2));
            re = zeros(size(Te,2),1);
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                re = re + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*bcInfo{i}{4}*(bcInfo{i}{3} - psi.sf(ig,:)*u(Te)) ...
                                    -dt*(1-alpha)*psi.sf(ig,:)'*bcInfo{i}{4}*(bcInfo{i}{3} - psi.sf(ig,:)*u_n(Te)))*det(J);
                Ke = Ke + gWts(ig)*(alpha*dt*(psi.sf(ig,:)'*bcInfo{i}{4}*psi.sf(ig,:)))*det(J);
            end 
            r(Te) = r(Te) + re;
            K(Te,Te) = K(Te, Te) + Ke;    
        end
        
    % Same procedure for radiative boundary      
    elseif bcInfo{i}{2} == 3
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 4]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            Ke = zeros(size(Te,2));
            re = zeros(size(Te,2),1);
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                re = re + gWts(ig)*(-alpha*dt*(psi.sf(ig,:)'*sig*bcInfo{i}{4}*(bcInfo{i}{3}^4 - (psi.sf(ig,:)*u(Te))^4)) ...
                                    -dt*(1-alpha)*psi.sf(ig,:)'*sig*bcInfo{i}{4}*(bcInfo{i}{3}^4 - (psi.sf(ig,:)*u_n(Te))^4))*det(J);
                Ke = Ke + gWts(ig)*(alpha*dt*(psi.sf(ig,:)'*4*sig*bcInfo{i}{4}*(psi.sf(ig,:)*u(Te))^3*psi.sf(ig,:)))*det(J);
            end
            r(Te) = r(Te) + re;
            K(Te,Te) = K(Te, Te) + Ke;    
        end
    end
end

end

