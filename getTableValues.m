function [var, dvar] = getTableValues(name, value)
%-----------------------------------------------------------------------------------
% Description: 
%
% Input Variables : 
%
% Output Variables : 
%-----------------------------------------------------------------------------------
dif = abs(name(:,1) - ones(length(name(:,1)),1)*value);

[d, aux] = min(dif);

if d~=0
    ds=name(aux,1)-value;
    if (ds>0 && aux~=1) || aux == length(name(:,1))
        x1=name(aux,1);
        x2=name(aux-1,1);
        y1=name(aux,2);
        y2=name(aux-1,2);
    else
        x1=name(aux,1);
        x2=name(aux+1,1);
        y1=name(aux,2);
        y2=name(aux+1,2);
    end
    var = y1-(y2-y1)/(x2-x1)*(x1-value);
else
    if aux == length(name(:,1))
        x1=name(aux,1);
        x2=name(aux-1,1);
        y1=name(aux,2);
        y2=name(aux-1,2);
    else
        x1=name(aux,1);
        x2=name(aux+1,1);
        y1=name(aux,2);
        y2=name(aux+1,2);
    end
    var = name(aux,2);
end

dvar=(y2-y1)/(x2-x1);

end