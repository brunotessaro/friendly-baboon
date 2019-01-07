function [psi, gWts, gPts] = bar2n()

gPts = [-1/sqrt(3); 1/sqrt(3)];
gWts = [1 1];
psi.sf = [1/2*(1-gPts) 1/2*(1+gPts)];
psi.dcsi = [-1/2*ones(size(gPts,1),1) 1/2*ones(size(gPts,1),1)];
psi.deta = [];

end