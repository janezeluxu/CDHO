function [] = writevtktau(taus,tauDC,taut_percent,tauadv_percent,taudiff_percent)

fileID = fopen('./testvtk/taus.txt','w');
formatSpec = '%2.16f \n';
[nrows,ncols] = size(taus);
value = zeros(nrows,1);
value(:,1) = taus;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

fileID = fopen('./testvtk/tauDC.txt','w');
formatSpec = '%2.16f \n';
[nrows,ncols] = size(tauDC);
value = zeros(nrows,1);
value(:,1) = tauDC;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

fileID = fopen('./testvtk/taut_percent.txt','w');
formatSpec = '%2.16f \n';
[nrows,ncols] = size(taut_percent);
value = zeros(nrows,1);
value(:,1) = taut_percent;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

fileID = fopen('./testvtk/tauadv_percent.txt','w');
formatSpec = '%2.16f \n';
[nrows,ncols] = size(tauadv_percent);
value = zeros(nrows,1);
value(:,1) = tauadv_percent;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

fileID = fopen('./testvtk/taudiff_percent.txt','w');
formatSpec = '%2.16f \n';
[nrows,ncols] = size(taudiff_percent);
value = zeros(nrows,1);
value(:,1) = taudiff_percent;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);
end