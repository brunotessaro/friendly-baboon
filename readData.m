function readData(pathName,fileName)

global Elements Sections Materials Nodes Lines Concentrated Data ElemData LineTraction

disp('Begin read data');
tic;

Data.PathName=pathName;
Data.FileName=fileName(1:(size(fileName,2)-4));

% open file
fid=fopen([pathName fileName]);

% read name of problem
aux = read_line(fid);
Data.name=aux;

% read type of simplification
aux = read_line(fid);
aux = sscanf(aux, '%s',2);   

if strcmpi(aux,'planestress')
    Materials.PlaneStress=true;
else
    Materials.PlaneStress=false;
end

% read analysis type
aux = read_line(fid);
aux = sscanf(aux, '%i');
Data.analysis_type=aux;

% read number of nodes, number of elements, number of materials, number of
% sections, number of concentrated loads, number of loaded sides, number of line results, 
% number of constraints
aux = read_line(fid);
aux = sscanf(aux, '%i %i %i %i %i %i %i %i');

Data.nnodes=aux(1);
Data.nelements=aux(2);
Data.nmaterials=aux(3);
Data.nsections=aux(4);
Data.nconc=aux(5);
Data.nlinetraction=aux(6);
Data.nlineout=aux(7);
Data.nconstraints=aux(8);  


% read nodal data
for i=1:Data.nnodes
    
    aux = read_line(fid);
    aux = sscanf(aux, '%i %g %g  %g %g %i %i %g %g');
    
    % coordinates
    Nodes(aux(1)).coord(1:2)=aux(2:3,1);
    % orientation
    Nodes(aux(1)).orientation(1:2,1)=aux(4:5,1);   
    % code for BC in the order: u1, u2
    Nodes(aux(1)).incid(1:2)=aux(6:7,1);
    % value of prescribed BC
    Nodes(aux(1)).F(1:2)=aux(8:9,1);
    
end

%
% element shape:
% --------------
%
% 1 - triangular
% 2 - rectangular
%
% element types:
% -------------
%
%    Element types
%
%     |           |        Number of nodes for interpolation
%Type | Element   |------------------------------------------------------
%     | Name      |Geometry  |          u
%
%  1 -  T3        (    3                3         ) 
%  2 -  T6        (    6                6         ) 
%  3 -  T10       (   10               10         ) 
%  4 -  T15       (   15               15         ) 
%  5 -  T21       (   21               21         ) 
%  6 -  QS4/QL4   (    4                4         ) 
%  7 -  QS8       (    8                8         ) 
%  8 -  QL9       (    9                9         ) 
%  9 -  QS12      (   12               12         ) 
% 10 -  QL16      (   16               16         ) 
% 11 -  QS16      (   16               16         ) 
% 12 -  QS20      (   20               20         ) 
% 13 -  QL25      (   25               25         ) 
% 14 -  QL36      (   36               36         ) 
%
% Element interpolation codes (Number of nodes)
%    
%   Triangular shapes
%
%     1 (3)          2 (6)          3 (10)         4 (15)         5 (21)
%
%     o              o              o              o              o            
%     | \            | \            | \            | \            | \          
%     |  \           |  \           |  \           o  o           o  o         
%     |   \          |   \          o   o          |   \          |   \        
%     |    \         |    \         |    \         |    \         o  o o       
%     |     \        o     o        |     \        o  o  o        |     \      
%     |      \       |      \       |      \       |      \       o  o o o     
%     |       \      |       \      o   o   o      |       \      |       \    
%     |        \     |        \     |        \     o  o  o  o     o  o o o o   
%     |         \    |         \    |         \    |         \    |         \  
%     o----------o   o-----o----o   o---o---o--o   o--o--o--o-o   o--o-o-o-o-o 
%

%
%   Rectangular shapes
%
%        1 (4)           2 (8)          3 (9)         4 (12)          5 (16)
%                                                                             
%     o----------o   o-----o----o   o-----o----o   o--o----o--o   o--o----o--o
%     |          |   |          |   |          |   |          |   |          |
%     |          |   |          |   |          |   |          |   |          |
%     |          |   |          |   |          |   o          o   o  o    o  o
%     |          |   |          |   |          |   |          |   |          |
%     |          |   o          o   o     o    o   |          |   |          |
%     |          |   |          |   |          |   |          |   |          |
%     |          |   |          |   |          |   o          o   o  o    o  o
%     |          |   |          |   |          |   |          |   |          |
%     |          |   |          |   |          |   |          |   |          |
%     o----------o   o-----o----o   o-----o----o   o--o----o--o   o--o----o--o

%        6 (16)         7 (20)         8 (25)         9 (36)      
%
%     o--o--o--o-o   o--o-o-o-o-o   o--o--o--o-o   o--o-o-o-o-o   
%     |          |   |          |   |          |   |          |  
%     o          o   o          o   o          o   o  o o o o o 
%     |          |   |          |   |          |   |          |  
%     |          |   o          o   |          |   o  o o o o o  
%     o          o   |          |   o          o   |          |  
%     |          |   o          o   |          |   o  o o o o o  
%     |          |   |          |   |          |   |          | 
%     o          o   o          o   o          o   o  o o o o o  
%     |          |   |          |   |          |   |          |   
%     o--o--o--o-o   o--o-o-o-o-o   o--o--o--o-o   o--o-o-o-o-o                                                                                             
%
%
% number of nodes for each interpolation 
%
% triangular
%         code   1 2  3  4  5
nnodes(1).value=[3 6 10 15 21];
%
% quadrilateral
%         code   1 2 3  4  5  6  7  8  9
nnodes(2).value=[4 8 9 12 16 16 20 25 36];
%
% Set Element Data for all element types
%
% T3
%
type=1;
ElemData(type).shape=1;     % geometry shape
ElemData(type).code=1;  % interpolation codes for geometry and u
ElemData(type).node.coord=[...% coordinates of all nodes in master element
    0 0     %  node 1
    1 0     %  node 2
    0 1];   %  node 3
%
% T6
%
type=2;
ElemData(type).shape=1;     
ElemData(type).code=2;  
ElemData(type).node.coord=[...
    0 0     
    1 0     
    0 1      
    1/2 0
    1/2 1/2
    0   1/2];   
% 
% T10
%
type=3;
ElemData(type).shape=1;     
ElemData(type).code=3;  
ElemData(type).node.coord=[...
    0 0 
    1 0 
    0 1 
    1/3 0 
    2/3 1/3 
    0 2/3 
    2/3 0 
    1/3 2/3 
    0 1/3 
    1/3 1/3];   
% 
% T15
%
type=4;
ElemData(type).shape=1;     
ElemData(type).code=4;  
ElemData(type).node.coord=[...
    0 0 
    1 0 
    0 1 
    1/4 0 
    3/4 1/4 
    0 3/4 
    1/2 0 
    1/2 1/2 
    0 1/2 
    3/4 0 
    1/4 3/4 
    0 1/4 
    1/4 1/4 
    1/2 1/4 
    1/4 1/2];   
% 
% T21
%
type=5;
ElemData(type).shape=1;     
ElemData(type).code=5; 
ElemData(type).node.coord=[...
    0 0 
    1 0 
    0 1 
    1/5 0 
    4/5 1/5 
    0 4/5 
    2/5 0 
    3/5 2/5 
    0 3/5 
    3/5 0 
    2/5 3/5 
    0 2/5 
    4/5 0 
    1/5 4/5 
    0 1/5 
    1/5 1/5 
    3/5 1/5 
    1/5 3/5 
    2/5 1/5 
    2/5 2/5 
    1/5 2/5];   
% 
% QS4/QL4
%
type=6;
ElemData(type).shape=2;       
ElemData(type).code=1;
ElemData(type).node.coord=[...
    -1 -1
     1 -1
     1  1
    -1  1]; 
% 
% QS8
%
type=7;
ElemData(type).shape=2;       
ElemData(type).code=2; 
ElemData(type).node.coord=[...
    -1 -1
     1 -1
     1  1
    -1  1
     0 -1
     1  0
     0  1
    -1  0]; 
%
% QL9
%
type=8;
ElemData(type).shape=2;       
ElemData(type).code=3;  
ElemData(type).node.coord=[...
    -1 -1
     1 -1
     1  1
    -1  1
     0 -1
     1  0
     0  1
    -1  0
     0  0]; 
%
% QS12
%
type=9;
ElemData(type).shape=2;       
ElemData(type).code=4; 
ElemData(type).node.coord=[...
    -1 -1 
    1 -1 
    1 1 
    -1 1 
    -(1/3) -1 
    1 -(1/3) 
    1/3 1 
    -1 1/3 
    1/3 -1 
    1 1/3 
    -(1/3) 1 
    -1 -(1/3)]; 
%
% QL16
%
type=10;
ElemData(type).shape=2;       
ElemData(type).code=5;  
ElemData(type).node.coord=[...
-1 -1 
1 -1 
1 1 
-1 1 
-(1/3) -1 
1 -(1/3) 
1/3 1 
-1 1/3 
1/3 -1 
1 1/3 
-(1/3) 1 
-1 -(1/3) 
-(1/3) -(1/3) 
1/3 -(1/3) 
1/3 1/3 
-(1/3) 1/3];
%
% QS16
%
type=11;
ElemData(type).shape=2;       
ElemData(type).code=6;
ElemData(type).node.coord=[...
-1 -1
1 -1
1 1
-1 1
-(1/2) -1
1 -(1/2)
1/2 1
-1 1/2
0 -1
1 0
0 1
-1 0
1/2 -1
1 1/2
-(1/2) 1
-1 -(1/2)];
%
% QS20
%
type=12;
ElemData(type).shape=2;       
ElemData(type).code=7; 
ElemData(type).node.coord=[...
-1 -1
1 -1
1 1
-1 1
-(3/5) -1
1 -(3/5)
3/5  1
-1 3/5
-(1/5) -1
1 -(1/5)
1/5 1
-1 1/5
1/5 -1
1 1/5
-(1/5) 1
-1 -(1/5)
3/5 -1
1, 3/5
-(3/5) 1
-1 -(3/5)];
%
% QL25
%
type=13;
ElemData(type).shape=2;       
ElemData(type).code=8; 
ElemData(type).node.coord=[...
-1 -1 
1 -1 
1 1 
-1 1 
-(1/2) -1 
1 -(1/2) 
1/2 1 
-1 1/2 
0 -1 
1 0 
0 1 
-1 0 
1/2 -1 
1 1/2 
-(1/2) 1 
-1 -(1/2)
-(1/2) -(1/2) 
1/2 -(1/2) 
1/2 1/2 
-(1/2) 1/2 
0 -(1/2) 
1/2 0 
0 1/2 
-(1/2) 0 
0 0];  
%
% QL36
%
type=14;
ElemData(type).shape=2;       
ElemData(type).code=9; 
ElemData(type).node.coord=[...
-1 -1 
1 -1 
1 1 
-1 1 
-(6/10) -1 
1 -(6/10) 
6/10 1 
-1 6/10 
-(2/10) -1 
1 -(2/10) 
2/10 1 
-1 2/10 
2/10 -1 
1 2/10 
-(2/10) 1 
-1 -(2/10) 
6/10 -1 
1 6/10 
-(6/10) 1 
-1 -(6/10) 
-(6/10) -(6/10) 
6/10 -(6/10) 
6/10 6/10 
-(6/10) 6/10 
-(2/10) -(6/10) 
6/10 -(2/10) 
2/10 6/10 
-(6/10) 2/10 
2/10 -(6/10) 
6/10 2/10 
-(2/10) 6/10 
-(6/10) -(2/10) 
-(2/10) -(2/10) 
2/10 -(2/10) 
2/10 2/10 
-(2/10) 2/10];

for itype=1:length(ElemData)
    ElemData(itype).nnodes=...
        nnodes(ElemData(itype).shape).value(ElemData(itype).code);
       
end

% read element data
for ielem=1:Data.nelements
    
    aux = read_line(fid);    
    aux1 = sscanf(aux, '%*s %s',1);    
    
    switch aux1
        case 'T3'
            type=1;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %g %g %g');
        case 'T6'
            type=2;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'T10'
            type=3;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'T15'
            type=4;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'T21'
            type=5;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case {'QL4','QS4'}
            type=6;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QS8'
            type=7;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QL9'
            type=8;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QS12'
            type=9;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QL16'
            type=10;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QS16'
            type=11;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QS20'
            type=12;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QL25'
            type=13;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
        case 'QL36'
            type=14;
            aux1 = sscanf(aux, '%i %*s %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %g %g %g');
            
    end
    
    % total number of nodes of the element
    ntnodes=ElemData(type).nnodes;
    % element type
    Elements(aux1(1)).type=type;
    % element total nodes
    Elements(aux1(1)).nodes(1:ntnodes)=aux1(2:(1+ntnodes),1);
    % integration rule for K, F_domain and M
    Elements(aux1(1)).intrule(1:3)=aux1(2+ntnodes:4+ntnodes,1);    
    % material number
    Elements(aux1(1)).mat(1)=aux1(5+ntnodes,1);
    % section number
    Elements(aux1(1)).sec(1)=aux1(6+ntnodes,1);    
    % uniform load value
    Elements(aux1(1)).unifload(1:2,1)=aux1(7+ntnodes:8+ntnodes,1);
    % uniform temperature value
    Elements(aux1(1)).uniftemp(1,1)=aux1(9+ntnodes,1);
    
end

% read material data
for i=1:Data.nmaterials
    
    aux = read_line(fid);
    aux = sscanf(aux, '%i %g %g %g %g');
    
    % E
    Materials(aux(1)).E(1)=aux(2,1);
    % niu
    Materials(aux(1)).niu(1)=aux(3,1);
    % alpha
    Materials(aux(1)).alpha(1)=aux(4,1);
    % ro
    Materials(aux(1)).ro(1)=aux(5,1);
    
end

% read section data
for i=1:Data.nsections
    
    aux = read_line(fid);
    aux = sscanf(aux, '%i %g');
    
    % h
    Sections(aux(1)).h(1)=aux(2,1);
    
    
end

% read concentrated loads data
for i=1:Data.nconc
    
    aux = read_line(fid);
    aux = sscanf(aux, '%i %g %g %g %g');
    
    % x
    Concentrated(aux(1)).x=aux(2:3);
    % value
    Concentrated(aux(1)).F(1:2,1)=aux(4:5);
    
end

% read lines applied tractions data
for i=1:Data.nlinetraction
    
    aux = read_line(fid);
    aux2 = sscanf(aux, '%i', 3);
        
    % number of nodes
    LineTraction(aux2(1)).nnodes=aux2(2);
    
    % integration rule
    LineTraction(aux2(1)).intrule=aux2(3);
    
    aux3 = sscanf(aux, '%i', aux2(2)+3);
    
    % list of nodes
    LineTraction(aux2(1)).nodes(1:aux2(2))=aux3(4:(3+aux2(2)));
    
    aux3 = sscanf(aux, '%g');    
    
    % prescribed value list of tractions
    LineTraction(aux2(1)).value(1:aux2(2),1:2)=...
        reshape(aux3(4+aux2(2):3+3*aux2(2)),[aux2(2) 2]);
end

% read line results data
for i=1:Data.nlineout
    
    aux = read_line(fid);
    aux = sscanf(aux, '%i %g %g %g %g');
    
    % (start point) x1,  y1
    Lines(aux(1)).noi=aux(2:3);
    % (end point) x2, y2
    Lines(aux(1)).nof=aux(4:5);
    
end

% read constraints data of the form
% d slave + c1 d master 1 + c2 d master 2 + ... - q1 = 0
for i=1:Data.nconstraints
    
    aux = read_line(fid);
    aux = sscanf(aux, '%g');
    
    % number of master combination dofs
    Constraint(aux(1)).master.n=aux(2);
    for imdof=1:aux(2)
        Constraint(aux(1)).master.data(i,1:3)=aux(2+3*(i-1)+1:2+3*(i-1)+3);
    end
    Constraint(aux(1)).value=aux(2+3*(aux(2)-1)+3+1);
end

% read membrane tensor moment data (only for linear stability analysis)
if Data.analysis_type==3
    
    aux = read_line(fid);
    aux = sscanf(aux, '%g %g %g');
    
    % membrane tensor (n11, n22, n12)
    % (the usual sign convention is used for n alpha beta)
    Sections.membrane_tensor=[aux(1) aux(3);aux(3) aux(2)];
end

tElapsed=toc;    
disp(['Finished read data in ' num2str(tElapsed) 's']);


% close file
fclose(fid);   
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aux=read_line(fid)

read=true;

while read
    aux = fgets(fid);

    if ~strcmp(aux(1,1),'#')
        read=false;
    end
end