function [u_xp] = solOnArbitraryPos(nodeInfo, elemInfo, u, xp)
%-----------------------------------------------------------------------------------
% Description: This function is desingned to obtain the solution of an arbitraty 
%              position point.
%               
% Input Variables : nodeInfo = struct containing node related info.
%                   elemInfo = struct containing element related info.
%                   u = problem solution (all steps and positions).
%
% Output Variables : u_xp = solution on specific position
%-----------------------------------------------------------------------------------

% Get mesh parameters
X = nodeInfo.X;
T = elemInfo.T;
elemType = elemInfo.elemType;
volElemIdx = elemInfo.volElemIdx;

for i=1:size(volElemIdx,1)
    
    % Definition of initial guess
    csi0 = 0.0001;
    
    % Retrieve element data
    iElem = volElemIdx(i);
    [psi, nElNds, dim] = elementCall2(elemType(iElem), csi0);
    Te = T(iElem,1:nElNds);
    Xe = X(Te,1:dim);
    
    % Calculate residual and tangent
    R = psi.sf*Xe - xp;
    dR = psi.dcsi*Xe;
    
    % Start of Newton loop
    while(norm(R) > 1e-8)
        
        % Calculate increment and uptdate cs0
        csi0 = csi0 - R/dR;
        
        % Get element info on new csi0
        [psi, ~, ~] = elementCall2(elemType(iElem), csi0);
        
        % Calculate residual and tangent
        R = psi.sf*Xe - xp;
        dR = psi.dcsi*Xe;
    end
    
    % Check if csi0 is inside element
    if (csi0 >= -1 && csi0 <= 1)
        
        % Get solution on element
        u_xp = psi.sf*u(Te,:);
        
        % Stop the search
        break
    end
    
    % If element is not found xp is outside of the domain
    if (i == size(volElemIdx,1))
        return 
    end
end

end