%------------------------------------------------------------------------------------------------------
% Matlab code for getting mesh and geometry info.
% Author: Bruno Aguirre Tessaro, IST - Lisboa.
% Contact: bruno.tessaro@gmail.com
% Date: 05/12/2018
%------------------------------------------------------------------------------------------------------
clear, clc, clear

%% Geometry input

% Volume, boundary conditions and tags explanation:
%
% tags.vol: tag representing volume.
% tags.boundary: tags representing boundary, one for each BC.
%                           
%   Physical Tag: tags.boundary{any}{011}
%   Type of BC: tags.boundary{any}{2} 
%   BC values: tags.boundary{any}{3:end} 
%
%   Types can be: {2} = 0 - imposed temp     ->  where -> {3} = imposed temp value 
%                 {2} = 1 - imposed flux     ->  where -> {3} = imposed flux value 
%                 {2} = 2 - convective flux  ->  where -> {3:4} = ambient temp and convctive coeff
%                 {2} = 3 - radiative flux   ->  where -> {3:4} = ambient temp and emissivety

%------------------------------------------------------------------------------------------------------

fileName = 'geoAndMesh/testCaseBathe.msh';
tags.vol = 100;
tags.boundary{1} = {10, 2, 0, 0.04};
tags.boundary{2} = {11, 3, 0, 1};
k = [0.01 0; 0 0.01];
Q = 0;
rho = 1;
Cp = 0.01;
u0 = 1498.1505;
A = 1;
eta = 1;

[nodeInfo, elemInfo, bcInfo] = getMeshInfo(fileName, tags);
nNds = size(nodeInfo.X,1);

%% Physical parameters

params.sig = 1.38064852e-23;
params.k = k;
params.Cp = Cp;
params.rho = rho;
params.A = A;
params.eta = eta;
params.Q = Q*ones(nNds,1);
u0 = u0*ones(size(nNds,1),1);

%% Time discretization

t0 = 0;
tf = 1;
dt = 0.1;
alpha = 1;
nStep = fix((tf-t0)/dt+1);
t = linspace(t0,tf,nStep);

params.dt = dt;
params.alpha = alpha;

%% Matrix initialization and BC/IC imposition

% Solution initialization
u = zeros(nNds,nStep);
c = zeros(nNds,nStep);
dphi = zeros(2*nNds,1);

% Initial condition imposition
u(:,1) = u0;
c(:,1) = 1;

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
    [rates] = getFieldRates(nodeInfo, elemInfo, bcInfo, u(:,n), c(:,n), params);
    
    % Initialize next newton step
    u(:,n+1) = u(:,n);
    c(:,n+1) = c(:,n);
    
    % Calculate residual and tangent matrix
    [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), c(:,n+1), c(:,n), rates, params);
    
    % Calculate aux parameters
    iter = 0;
    ri = norm(r(nodeInfo.free));
    
    % Start newton loop
    while(norm(r(nodeInfo.free))/ri > 1e-6)
        
        % Iteration counter
        iter=iter+1;
            
        % Calculate solution on the next iterative step
        dphi(nodeInfo.free) = - K(nodeInfo.free,nodeInfo.free)\r(nodeInfo.free);
        
        % Mount solution vectors
        u(:,n+1) =  u(:,n+1) + dphi(1:nNds);
        c(:,n+1) =  c(:,n+1) + dphi(nNds+1:2*nNds);
        
        % Calculate residual and tangent matrix
        [K, r] = elementSubRoutine(nodeInfo, elemInfo, bcInfo, u(:,n+1), u(:,n), c(:,n+1), c(:,n), rates, params);
        
        % Display iteration info
        disp(['Time = ' num2str(dt*n)  '  Iter = ' num2str(iter) '    Res = ' num2str(norm(r(nodeInfo.free))) '    ResAdm = ' num2str(norm(r(nodeInfo.free))/norm(ri))]);     
        
    end
end

%% Writting gmsh files for vizualization
writeGmsh(fileName, nodeInfo, u, t)
writeGmsh(fileName, nodeInfo, c, t)


