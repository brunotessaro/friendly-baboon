function [absc, wght] = gaussIntrgParams(ord)
%% Generates the abscissa and weights for a Gauss-Legendre quadrature. 
if (ord==1)
  absc = 0;
  wght = 2;
  return
end
vect = (1:ord-1)./sqrt(4*((1:ord-1)).^2-1);
[w,x] = eig(diag(vect,-1)+diag(vect,1));
absc = diag(x)';
wght = 2*(w(1,:)').^2;
