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
nNds = size(X,1);

% Get physical parameters
k = params.k;
Q = params.Q;
sig = params.sig;
rho = params.rho;
Cp = params.Cp;
dt = params.dt;
alpha = params.alpha;

% Initialize global tangents and resdiuals
r_u = zeros(nNds,1);
r_c = zeros(nNds,1);
K_uu = zeros(nNds);
K_uc = zeros(nNds);
K_cu = zeros(nNds);
K_cc = zeros(nNds);


%% Assembly of volume elements contributions

% Loop trough volume elements
for i=1:size(volElemIdx,1)
    
    % Retrieve element data and initialize elemental matrices
    iElem = volElemIdx(i);
    [psi, gWts, nElNds, dim] = elementCall(elemType(iElem), [2 2]);
    
    Te = T(iElem,1:nElNds);
    Xe = X(Te,:);
    re_u = zeros(size(Te,2),1);
    re_c = zeros(size(Te,2),1);
    Ke_uu = zeros(size(Te,2));
    Ke_cu = zeros(size(Te,2));
    Ke_uc = zeros(size(Te,2));
    Ke_cc = zeros(size(Te,2));
    
    % Loop trough element gauss points
    for ig=1:size(gWts,2)
            J = jacobCalc(elemType(iElem), ig, Xe, psi);
            
            % Assemble B matrix depending on dimension.
            if dim == 1
                B = J\psi.dcsi(ig,:); 
            elseif dim == 2
                B = J\[psi.dcsi(ig,:) ; psi.deta(ig,:)];
            end
            
%%%---------------------- Arrange later ----------------------%%          
            rho = 1;
            rho_n = 1;
            Cp = 1;
            Cp_n = 1;
            k = 1;
            k_n = 1;
%%%-----------------------------------------------------------%%          

    
            M = psi.sf(ig,:)'*psi.sf(ig,:);
            Kappa = B'*B;
            F_Q = M*Q(Te);
            
            re_u = re_u + gWts(ig)*(M*rho*Cp*(u(Te) - u_n(Te)) ...
                                + alpha*dt*(Kappa*k*u(Te) - F_Q) ...
                                + dt*(1-alpha)*(rho*Cp)/(rho_n*Cp_n)*(Kappa*k_n*u_n(Te) - F_Q))*det(J);     
            
            %Ke = %Ke + gWts(ig)*(M + alpha*dt*Kappa)*det(J);
            
    end 
    
    % Assembly of residual and tangent matrix
     r_u(Te) = r_u(Te) + re_u;
     r_c(Te) = r_c(Te) + re_c;
     K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;
     K_uc(Te,Te) = K_uc(Te, Te) + Ke_uc;
     K_cu(Te,Te) = K_cu(Te, Te) + Ke_cu;
     K_cc(Te,Te) = K_cc(Te, Te) + Ke_cc;

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
            re_u = zeros(size(Te,2),1);
            
            % Loop trough element gauss points
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
              
%%%---------------------- Arrange later ----------------------%%          
            q = bcInfo{i}{3};
            q_n = bcInfo{i}{3};
%%%-----------------------------------------------------------%%   
               

                re_u = re_u + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*q ...
                                        -dt*(1-alpha)*(rho*Cp)/(rho_n*Cp_n)*psi.sf(ig,:)'*q_n)*det(J);
            
            end
            
            % Assembly of residual and tangent matrix
            r_u(Te) = r_u(Te) + re_u;            
        end
        
        % Same procedure for convective boundary
    elseif bcInfo{i}{2} == 2
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 1]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
            re_u = zeros(size(Te,2),1);
            Ke_uu = zeros(size(Te,2));
            Ke_uc = zeros(size(Te,2));
            
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
%%%---------------------- Arrange later ----------------------%%          
            h = bcInfo{i}{4};
            h_n = bcInfo{i}{4};
            u_a = bcInfo{i}{3};
%%%-----------------------------------------------------------%%  
                
                
                re_u = re_u + gWts(ig)*(-alpha*dt*psi.sf(ig,:)'*h*(u_a - psi.sf(ig,:)*u(Te)) ...
                                        -dt*(1-alpha)*(rho*Cp)/(rho_n*Cp_n)*psi.sf(ig,:)'*h_n*(u_a - psi.sf(ig,:)*u_n(Te)))*det(J);
                
                
                %Ke = %Ke + gWts(ig)*(alpha*dt*(psi.sf(ig,:)'*bcInfo{i}{4}*psi.sf(ig,:)))*det(J);
            end
            
            r_u(Te) = r_u(Te) + re_u;
            K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;
            K_uc(Te,Te) = K_uc(Te, Te) + Ke_uc;            
        end
        
        % Same procedure for radiative boundary
    elseif bcInfo{i}{2} == 3
        for j=1:size(bcInfo{i}{1},1)
            iElem = bcInfo{i}{1}(j);
            [psi, gWts, nElNds, ~] = elementCall(elemType(iElem), [1 4]);
            Te = T(iElem,1:nElNds);
            Xe = X(Te,:);
            
%%%---------------------- Arrange later ----------------------%%          
            eps = bcInfo{i}{4};
            eps_n = bcInfo{i}{4};
            u_a = bcInfo{i}{3};
%%%-----------------------------------------------------------%%  
            
            re_u = zeros(size(Te,2),1);
            Ke_uu = zeros(size(Te,2));
            Ke_uc = zeros(size(Te,2));
            
            for ig=1:size(gWts,2)
                J = jacobCalc(elemType(iElem), ig, Xe, psi);
                
                re_u = re_u + gWts(ig)*(-alpha*dt*(psi.sf(ig,:)'*sig*eps*(u_a^4 - (psi.sf(ig,:)*u(Te))^4)) ...
                                        -dt*(1-alpha)*(rho*Cp)/(rho_n*Cp_n)*psi.sf(ig,:)'*sig*eps_n*(u_a^4 - (psi.sf(ig,:)*u_n(Te))^4))*det(J);
                
                
                %Ke = %Ke + gWts(ig)*(alpha*dt*(psi.sf(ig,:)'*4*sig*bcInfo{i}{4}*(psi.sf(ig,:)*u(Te))^3*psi.sf(ig,:)))*det(J);
            end
            
            r_u(Te) = r_u(Te) + re_u;
            K_uu(Te,Te) = K_uu(Te, Te) + Ke_uu;
            K_uc(Te,Te) = K_uc(Te, Te) + Ke_uc;
        end
    end
end

end

