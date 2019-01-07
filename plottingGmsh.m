function [mapp] = plottingGmsh(T, X, u)

Taux = zeros(size(T));
for i=1:size(T,1)
    n = size(T(i,:),2);
    Taux(i,:) = [T(i,1) T(i,3:end) T(i,2)];       
end

mapp = zeros(size(X,2),1);
k = 0;
for i=1:size(Taux,1)
    for j=1:size(Taux,2)
        if(not(ismember(Taux(i,j),mapp)))
           k = k +1;
           mapp(k) = Taux(i,j);     
        end
    end
end
    
%plot(X(mapp,1),u(mapp))

end