function [] = plot1D(u, X, t, limy, limx)
%-----------------------------------------------------------------------------------
% Description: This function creates the residual and tangent matrix of the problem. 
%               
% Input Variables : u = solution matrix.
%                   X = space vector.
%                   t = time vector.
%
% Output Variables : -
% 
%-----------------------------------------------------------------------------------

% Get auxiliary vectors
nNds = size(X,1);
nStep = size(t,2);

% Arranging matrices for plotting (in 1D gmsh puts the last node in the second position)
u(nNds+1,:) = u(2,:);
u(2,:) = [];

X = X(:,1);
X(nNds+1,:) = X(2,:);
X(2,:) = [];

% Plot in space for each time step
figure(1)
for n=1:1:nStep
    plot(X,u(:,n))
    ylim(limy)
    xlim(limx)
    title(['Time: ', num2str(t(n)), 'min'])
    legend('Numerical Solution')
    pause(0.2)
end

end

