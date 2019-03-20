%-------------------------------------------------------------------------------------------------------
% Matlab code for solving the Bai's equations.
% Author: Bruno Aguirre Tessaro, IST - Lisboa.
% Contact: bruno.tessaro@tecnico.ulisboa.com
%-------------------------------------------------------------------------------------------------------
clear, clc, clear

%% Data reading and assignments

% Read data from input file
[matParams, meshParams, timeParams] = readData('inputs/NC-SLC03.inp');

% Get mesh information and number of nodes
[nodeInfo, elemInfo, bcInfo] = getMeshInfo(meshParams.fileName, meshParams.tags);
nNds = size(nodeInfo.X,1);

% Struct assigning and data fixing (here dt and alpha are introduced in params to save space)
params = matParams;
params.dt = timeParams.dt;
params.alpha = timeParams.alpha;
params.GFRP = load('GFRP.mat');

%% Time discretization

% Get number of steps
nStep = fix((timeParams.tf-timeParams.t0)/timeParams.dt + 1);

% Get time vector
t = linspace(timeParams.t0,timeParams.tf,nStep);

%% Calculation of time dependet variables

% Fire temperature
u_inf = timeParams.u0 + 345*log10(8*t/60+1);

%% Matrix initialization and BC/IC imposition

% Solution initialization
u = zeros(nNds,nStep);

% Initial condition imposition
u(:,1) = timeParams.u0*ones(nNds,1);

% Imposing temperature BC 
for i=1:size(bcInfo,2)
    for j=1:size(bcInfo{i}{1},1)
        if bcInfo{i}{2} == 0
            u(nonzeros(elemInfo.T(bcInfo{i}{1}(j),:)),:) = bcInfo{i}{3};
        end
    end
end

%% Start the time loop
for n=1:nStep-1
    
    % Calculate rates at first time step
    % Obs: Since rates are calculated in n and the n-1 is required for the calculation of Cp in the
    % material point one has to countermeasure when n == 1
    
    [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u(:,n), u_inf(n), params);
    
    % Initialize next newton step
    u(:,n+1) = u(:,n);
    
    % Calculate residual and tangent matrix
    [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), rates, u_inf(n+1), params);
    
    % Calculate aux parameters
    iter = 0;
    ri = norm(r(nodeInfo.free));
    
    % Display first iteration info
    disp(['Time = ' num2str(timeParams.dt*n)  '  Iter = ' num2str(iter) '    Res = ' num2str(norm(r(nodeInfo.free))) '    ResAdm = ' num2str(norm(r(nodeInfo.free))/ri)]);
    
    % Start newton loop
    while(norm(r(nodeInfo.free))/ri > 1e-6)
        
        % Iteration counter
        iter=iter+1;
        dphi = zeros(nNds,1);

        % Calculate solution on the next iterative step
        dphi(nodeInfo.free) = - K(nodeInfo.free,nodeInfo.free)\r(nodeInfo.free);
        
        % Mount solution vectors
        u(:,n+1) =  u(:,n+1) + dphi(1:nNds);
        
        % Calculate residual and tangent matrix
        [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), rates, u_inf(n+1), params);
        
        % Display iteration info
        disp(['Time = ' num2str(timeParams.dt*n)  '  Iter = ' num2str(iter) '    Res = ' num2str(norm(r(nodeInfo.free))) '    ResAdm = ' num2str(norm(r(nodeInfo.free))/norm(ri))]);     
        
    end
end

%% Writting gmsh files for vizualization
writeGmsh(meshParams.fileName, nodeInfo, u, t)


idx1 = find(abs(nodeInfo.X - 0.0163*0.25) < 0.0001);
idx2 = find(abs(nodeInfo.X - 0.0163*0.5) < 0.0001);
idx3 = length(nodeInfo.X);

idx = [idx1, idx2, idx3];

%Plotting in 1D
figure(1)
plot(t/60,u(idx1,:),t/60,u(idx2,:),t/60,u(idx3,:))
grid on

% figure(2)
% hold on
% grid on
% plot(t/60,c(end,:), '-.m')

% figure(2)
% hold on
% plot(t/60,c(2,:))

% plot1D(u, nodeInfo.X, t, [0 1000], [0 0.016])




