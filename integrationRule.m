function [absc, wght]=integrationRule(shape,intrule)

switch shape
    case 1
        [absc, wght] = Triangular(intrule);
        
    case 2
        ngauss=intrule;
        [absc1D, wght1D] = gaussIntrgParams(ngauss);
        
        for i=1:ngauss
            for j=1:ngauss
                absc(1:2,(i-1)*ngauss+j)=[absc1D(i); absc1D(j)];
                wght((i-1)*ngauss+j)=wght1D(i)*wght1D(j);
            end
        end
        
    case 3
        ngauss=intrule;
        [absc, wght] = gaussIntrgParams(ngauss);
        wght =  wght';
end

end