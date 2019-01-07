function [psi, gWts, gPts] = quad4n()

% gPts = [-sqrt(3/5); 0; sqrt(3/5)];
% csi1D = [-sqrt(3/5); 0; sqrt(3/5)];
% gWts1D = [5/9 8/9 5/9];
% 
% 
% 
% for i=1:3
%     for j=1:3
% csi(3*(i-1)+j) = csi1D(i);
% eta(3*(i-1)+j) = csi1D(j);
%     end
% end
% 
% gPts = [csi' eta'];
% csi = csi';
% eta = eta';
% 
% tensorprod=gWts1D'*gWts1D;
% 
% for i=1:3
%     for j=1:3
% gWts(3*(i-1)+j) =tensorprod(i,j);
% 
%     end
% end
% clear tensorprod;


csi = [-1/sqrt(3); 1/sqrt(3); 1/sqrt(3); -1/sqrt(3)];
eta = [-1/sqrt(3); -1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];

gPts = [csi eta];

gWts = [1 1 1 1];


% 4 node Serendipian/Lagrangian
                
%       4--------3
%       |        |
%       |        |
%       |        |
%       1--------2

% Shape functions
psi.sf = [(1/4).*((-1)+eta).*((-1)+csi), ...
          (-1/4).*((-1)+eta).*(1+csi), ...
          (1/4).*(1+eta).*(1+csi), ...
          (-1/4).*(1+eta).*((-1)+csi)];

% Derivatives wrt csi
psi.dcsi = [(1/4).*((-1)+eta), ...
           (1/4).*(1+(-1).*eta), ...
           (1/4).*(1+eta), ...
           (1/4).*((-1)+(-1).*eta)];

% Derivatives wrt eta
psi.deta= [(1/4).*((-1)+csi), ...
           (1/4).*((-1)+(-1).*csi), ...
           (1/4).*(1+csi), ... 
           (1/4).*(1+(-1).*csi)];

end