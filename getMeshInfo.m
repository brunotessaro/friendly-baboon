%-----------------------------------------------------------------------------------
% Description: This function loads the gmsh originated mesh and gathers parameters
%              relevant to the problem resolution.
%               
% Input Variables : fileName = mesh file name to be inputed.
%                   tags = struct containing physical tags from boundary and volume.
%
% Output Variables : nodeInfo = struct containing nodes information.
%                    nodeInfo = struct containing elements information.
%                    bcInfo = struct containing boundary condition information.
%-----------------------------------------------------------------------------------
function [nodeInfo, elemInfo, bcInfo] = getMeshInfo(fileName, tags)
% Load mesh                 
m = load_gmsh2(fileName,-1);

% Get position matrix and connectivety
X = m.POS;
T = m.ELE_NODES;

% Get element types
elemType = m.ELE_INFOS(:,2);

% Get ELEMENTAL indexes corresponding to volume
volElemIdx = find(m.ELE_TAGS(:,1) == tags.vol(1));

% Get integration rule for volume elements
volIntRule = tags.vol(2);

% Get boundary ELEMENTAL indexes 
for i=1:size(tags.boundary,2)
    bcElemIdx{i} = find(m.ELE_TAGS(:,1) == tags.boundary{i}{1});
end

% Get boundary NODE indexes and restrained nodes
restr = [];
for i=1:size(tags.boundary,2)
    bcIdx{i} = unique(nonzeros(T(bcElemIdx{i},:)), 'stable'); %#ok<*AGROW,*FNDSB>
    if tags.boundary{i}{2} == 0
        restr =  [restr ; bcIdx{i}];
    end
end

% Get NODE indexes corresponding to volume
volIdx = unique(nonzeros(T(volElemIdx,:)),'stable');

% Get free Nodes
free = volIdx;
for i=1:size(restr,1)
    free(find(free == restr(i))) = [];
end

% Get ELEMENTAL boundary information
for i=1:size(tags.boundary,2)
    tags.boundary{i}{6} = bcElemIdx{i};
end
bcInfo = tags.boundary;

% nodeInfo struct assignment
nodeInfo.X = X;
nodeInfo.free = free;
nodeInfo.restr = restr;
nodeInfo.volIdx = volIdx; 

% elemInfo struct assignment
elemInfo.T = T;
elemInfo.elemType = elemType;
elemInfo.volElemIdx = volElemIdx;
elemInfo.volIntRule = volIntRule;

end