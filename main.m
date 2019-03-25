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
params.Q = params.Q*ones(nNds,1);
params.dt = timeParams.dt;
params.alpha = timeParams.alpha;

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
c = zeros(nNds,nStep);
dphi = zeros(2*nNds,1);

% Initial condition imposition
u(:,1) = timeParams.u0*ones(nNds,1);
c(:,1) = timeParams.c0*ones(nNds,1);

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
    if n == 1
        [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u(:,n), u(:,n), c(:,n) ,c(:,n), u_inf(n), params);
    else
        [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u(:,n), u(:,n-1), c(:,n), c(:,n-1), u_inf(n+1), params);
    end
    % Initialize next newton step
    u(:,n+1) = u(:,n);
    c(:,n+1) = c(:,n);
    
    % Calculate residual and tangent matrix
    [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), c(:,n+1), c(:,n), rates, u_inf(n+1), params);
    
    % Calculate aux parameters
    iter = 0;
    ri = norm(r(nodeInfo.free));
    
    % Display first iteration info
    disp(['Time = ' num2str(timeParams.dt*n)  '  Iter = ' num2str(iter) '    Res = ' num2str(norm(r(nodeInfo.free))) '    ResAdm = ' num2str(norm(r(nodeInfo.free))/ri)]);
    
    % Start newton loop
    while(norm(r(nodeInfo.free))/ri > 1e-8)
        
        % Iteration counter
        iter=iter+1;
        dphi = zeros(2*nNds,1);

        % Calculate solution on the next iterative step
        dphi(nodeInfo.free) = - K(nodeInfo.free,nodeInfo.free)\r(nodeInfo.free);
        
        % Mount solution vectors
        u(:,n+1) =  u(:,n+1) + dphi(1:nNds);
        c(:,n+1) =  c(:,n+1) + dphi(nNds+1:2*nNds);
        
        % Since 'c' cannot be greater than 1, once it reaches this values the DoF is removed from the
        % system of equations (altering nodeInfo.free) and it is fixed to 1.
        aux = find(c(:,n+1)>1);
        if(~isempty(aux))       
            c(aux,n+1) = 1;
            for i=1:length(aux)
                nodeInfo.free(find(nodeInfo.free == aux(i)+nNds)) = [];
            end
        end
        
        % Calculate residual and tangent matrix
        [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), c(:,n+1), c(:,n), rates, u_inf(n+1), params);
        
        % Display iteration info
        disp(['Time = ' num2str(timeParams.dt*n)  '  Iter = ' num2str(iter) '    Res = ' num2str(norm(r(nodeInfo.free))) '    ResAdm = ' num2str(norm(r(nodeInfo.free))/norm(ri))]);     
        
    end
end

%% Writting gmsh files for vizualization
writeGmsh(meshParams.fileName, nodeInfo, u, t)
writeGmsh(meshParams.fileName, nodeInfo, c, t)

%% Get solution on arbitraty positions
xp1 = 4.2/1000;
u_xp1 = solOnArbitraryPos(nodeInfo, elemInfo, u, xp1);

xp2 = 0.0162;
u_xp2 = solOnArbitraryPos(nodeInfo, elemInfo, u, xp2);

xp3 = 0.0162;
u_xp2 = solOnArbitraryPos(nodeInfo, elemInfo, u, xp2);

%% Plotting in 1D
% 
figure(1)
plot(t/60,u_xp1)
legend('d = L/4 (mm)', 'd = L/2 (mm)', 'd = L (mm)', 'Location', 'northwest')
ylim([0 800])
grid on
hold on

% WF_dt15_n65 = u_xp1;
% WF_dt15_n65_x_end = u_xp2;
% save('../results/WF_dt15_n65.mat', 'WF_dt15_n65')
% save('../results/WF_dt15_n65_x_end.mat', 'WF_dt15_n65_x_end')





