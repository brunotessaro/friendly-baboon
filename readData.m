function [matParams, meshParams, timeParams] = readData(fileName)
%-----------------------------------------------------------------------------------
% Description: This function reads the data written on an input file. 
%               
% Input Variables : fileName = name of input (.inp) file.
%
% Output Variables : matParams = struct containing material parameters
%                    meshParams = struct containing mesh parameters
%                    timeParams = struct containing time parameters
%-----------------------------------------------------------------------------------
%% Read data and assign to corresponding struct

disp('Begin read data');
tic;

% open file
fid=fopen(fileName);

% read name of problem
aux = read_line(fid);
params.name=aux;

% read mesh name
aux = sscanf(read_line(fid), '%s');
meshParams.fileName = aux;

% read volume physical tag
aux = sscanf(read_line(fid), '%i %i');
meshParams.tags.vol(1) = aux(1);
meshParams.tags.vol(2) = aux(2);

% read boundary conditions
aux = sscanf(read_line(fid), '%i');
nPhysiTags = aux;
for i=1:nPhysiTags
    aux = sscanf(read_line(fid), '%i %i %g %g %i');
    for j=1:5
        meshParams.tags.boundary{i}{j} = aux(j);
    end
end

% read time parameters
aux = sscanf(read_line(fid), '%g %g %g %i');
timeParams.t0 = aux(1);
timeParams.tf = aux(2);
timeParams.dt = aux(3);
timeParams.alpha = aux(4);

% read initial conditions
aux = sscanf(read_line(fid), '%g %g');
timeParams.u0 = aux(1);
timeParams.c0 = aux(2);

% read general physical parameters
aux = sscanf(read_line(fid), '%g %g %g %g %g %g %g');
matParams.E = aux(1);
matParams.A = aux(2);
matParams.eta = aux(3);
matParams.R = aux(4);
matParams.Cp_d = aux(5);
matParams.sig = aux(6);
matParams.Q = aux(7);

% read before and after decomposition physical parameters
aux = sscanf(read_line(fid), '%g %g %g %g %g %g');
matParams.rho_b = aux(1);
matParams.rho_a = aux(2);
matParams.Cp_b = aux(3);
matParams.Cp_a = aux(4);
matParams.k_b = aux(5);
matParams.k_a = aux(6);

% read boundary physical parameters
aux = sscanf(read_line(fid), '%g %g %g %g %g %g %g');
matParams.k_g = aux(1);
matParams.Pr = aux(2);
matParams.g = aux(3);
matParams.beta = aux(4);
matParams.nu = aux(5);
matParams.eps_1 = aux(6);
matParams.eps_2 = aux(7);

fclose(fid);
disp('End read data');
toc;
% -------------------------------------------------------------------------------------------------------

function aux=read_line(fid)
%% Auxiliary function for line reading

read=true;

while read
    aux = fgets(fid);

    if ~strcmp(aux(1,1),'#')
        read=false;
    end
end
