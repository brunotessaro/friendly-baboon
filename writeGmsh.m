function [] = writeGmsh(fileName, nodeInfo, u, t)


disp('Begin write gmsh file');
tic;

% Copy file
copyfile(fileName, 'results.msh')  

% Open file
fid = fopen('results.msh','a+');

for i=1:length(t)-1
    % write displacements
    fprintf(fid,'$NodeData\n');
    
    % number of string tags
    fprintf(fid,'%13.0f\n',1);
    
    % string tags
    fprintf(fid,'"Temperature"\n');
    
    % number of real tags
    fprintf(fid,'%13.0f\n',1);
    % real tags
    fprintf(fid,'%27.15f\n',1.0); % time
    
    % number of integer tags
    fprintf(fid,'%13.0f\n',3);
    % integer tags
    fprintf(fid,'%13.0f\n',i); % time step index
    fprintf(fid,'%13.0f\n',1); % number of field components
    fprintf(fid,'%13.0f\n',size(nodeInfo.X,1)); % number of entities in the view
    
    for inode=1:size(nodeInfo.X,1)
        fprintf(fid,'%13.0f ',inode);
        fprintf(fid,'%25.15f ',u(inode,i));
        fprintf(fid,'\n');
    end
    fprintf(fid,'$EndNodeData\n');
end

% close file
fclose(fid);

tElapsed=toc;
disp(['Finished write gmsh file in ' num2str(tElapsed) 's';]);



