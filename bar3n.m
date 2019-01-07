function [psi, gWts, gPts] = bar3n()

gPts = [-sqrt(3/5); 0; sqrt(3/5)];
csi = [-sqrt(3/5); 0; sqrt(3/5)];
gWts = [5/9 8/9 5/9];

psi.sf = [(1/2).*((-1)+csi).*csi, (1/2).*csi.*(1+csi), 1+(-1).*csi.^2];
psi.dcsi = [(-1/2)+csi, (1/2)+csi, (-2).*csi];
psi.deta = [];

end