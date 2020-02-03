function [] = writeMesh(meshData,vertexData,file)
%fileID = fopen('./testvtk/cases/data/elements.txt','w');
filename = strcat('./',file,'.m');
fileID = fopen(filename,'w');
string1 = strcat('function [x,y,ien] = ',file,'() \n');
fprintf(fileID, string1);
fprintf(fileID, 'ien = [');
formatSpec = '%d %d %d\n';
[nrows,ncols] = size(meshData);
for row = 1:nrows
    fprintf(fileID,formatSpec,meshData{row,2});
end
fprintf(fileID, '];\n');
fclose(fileID);

fileID = fopen(filename,'a');
fprintf(fileID, 'x = [');
formatSpec = '%2.16f  \n';
[nrows,ncols] = size(vertexData);
for row = 1:nrows
    fprintf(fileID,formatSpec,vertexData{row,2}(1));
end
fprintf(fileID, '];\n');

fprintf(fileID, 'y = [');
formatSpec = '%2.16f  \n';
[nrows,ncols] = size(vertexData);
for row = 1:nrows
    fprintf(fileID,formatSpec,vertexData{row,2}(2));
end
fprintf(fileID, '];\n');
fprintf(fileID, 'end \n');
fclose(fileID);