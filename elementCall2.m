function [psi, gWts, nElNds, dim] = elementCall2(elmType)

% 2-noded bar element
if elmType == 1
    [psi, gWts, gPts] = bar2n();
    nElNds = size(psi.sf,2);
    dim = 1;

% 3-noded bar element
elseif elmType == 8
    [psi, gWts, gPts] = bar3n();
    nElNds = size(psi.sf,2);
    dim = 1;

% 4-noded quad element    
elseif elmType == 3
    [psi, gWts, gPts] = quad4n();
    nElNds = size(psi.sf,2);
    dim = 2;

% 1-noded point element    
elseif elmType == 15
    psi.sf = 1;
    gWts = 1;
    gPts = 1;
    nElNds = 1;
    dim = 0;

end