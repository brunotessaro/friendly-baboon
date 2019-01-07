function [J] = jacobCalc(elemType, ig, Xe, psi)
%-----------------------------------------------------------------------------------
% Description: This function calculates de jacobian matrix for each element. 
%               
% Input Variables : elemType = type of element. 
%                   ig = integration point index.
%                   Xe = element node positions.
%                   psi = shape function information.
%
% Output Variables : J = jacobian matrix.
%-----------------------------------------------------------------------------------
%% Calls jacobian matrix depending on element type 

% Point elements
if(elemType == 15)
     J = 1;
     
% Bar elements
elseif(elemType == 1 || elemType == 8)
     J = sqrt((psi.dcsi(ig,:)*Xe(:,1))^2  + (psi.dcsi(ig,:)*Xe(:,2))^2);

% Quadrilateral elements
elseif(elemType == 3 || elemType == 10)
      J = [psi.dcsi(ig,:)*Xe(:,1)  psi.dcsi(ig,:)*Xe(:,2); 
           psi.deta(ig,:)*Xe(:,1)  psi.deta(ig,:)*Xe(:,2)];
      
end

end