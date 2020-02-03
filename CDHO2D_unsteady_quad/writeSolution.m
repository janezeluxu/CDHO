function [] = writeSolution(u,filename)
fileID = fopen(filename,'w');
formatSpec = '%2.16f \n';
fprintf(fileID,formatSpec,u);
fclose(fileID);
end