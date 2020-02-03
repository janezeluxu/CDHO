function [Grid_size,nodeSize,totalDOF,pOrder,meshIEN,solution_ien,order_list,...
    edge_usage,edgeNode, vertexData,BCarray,BCele] = readinMesh(fileName)
%fileName = "./mesh/test_mesh.txt";
fileID = fopen(fileName,'r');
tline = fgetl(fileID);
i = 1;
while ischar(tline)
    C = textscan(tline,'%s');
    c= C{1};
    
    if ( strcmp(c{1},'Grid_size')==1)
        Grid_size = str2num(c{2});
    elseif ( strcmp(c{1},'Node_size')==1)
        nodeSize = str2num(c{2});
    elseif ( strcmp(c{1},'Total_DOF')==1)
        totalDOF = str2num(c{2});
    elseif ( strcmp(c{1},'maxDOF')==1)
        maxDOF = str2num(c{2});
    elseif ( strcmp(c{1},'Solution_Order')==1)
        npstart = i;
    elseif ( strcmp(c{1},'bcElement')==1)
        numBC = str2num(c{2});
    elseif ( strcmp(c{1},'edge_Size')==1)
        edgeNum = str2num(c{2});
    elseif ( strcmp(c{1},'edgeUsage')==1)
        nedgeUstart = i;
        %pOrder = str2num(order);
    end
    if ( strcmp(c{1},'meshIEN')==1)
        nIENStart = i;
    end
    
    if ( strcmp(c{1},'solutionIEN')==1)
        nsolutionStart = i;
    end
    
    if ( strcmp(c{1},'solutionOrder')==1)
        norderStart = i;
    end
    
    if ( strcmp(c{1},'vertexCoord')==1)
        nVStart = i;
    end
    
    if ( strcmp(c{1},'bc_tag')==1)
        nVcStart = i;
    end

    if ( strcmp(c{1},'EdgeNode')==1)
        nEdgeNodeStart = i;
    end
    
    if ( strcmp(c{1},'bcElement')==1)
        bceStart = i;
    end
    
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);

%celldisp(Array)
pOrder = str2num(Array{npstart+1});

meshIEN = zeros(Grid_size,4);
for nFace = 1:Grid_size
    meshIEN(nFace,:) = str2num(Array{nIENStart+nFace});
end

solution_ien = zeros(Grid_size,maxDOF);
for nFace = 1:Grid_size
    solution_ien(nFace,:) = str2num(Array{nsolutionStart+nFace});
end

order_list = zeros(Grid_size,5);
for nFace = 1:Grid_size
    order_list(nFace,:) = str2num(Array{norderStart+nFace});
end

edge_usage = zeros(Grid_size,3);
for nFace = 1:Grid_size
    edge_usage(nFace,:) = str2num(Array{nedgeUstart+nFace});
end

vertexData = zeros(nodeSize,2);
for nV = 1:nodeSize
    vertexData(nV,:) = str2num(Array{nVStart+nV});
end

%nodeSizet = 83;
BCarray = zeros(totalDOF,1);
for idV = 1:totalDOF
    BCarray(idV,:) = str2num(Array{nVcStart+idV});
end

maxOrder = max(pOrder);
edgeNode = zeros(edgeNum,maxOrder+1);
for idV = 1:edgeNum
    edgeNode(idV,:) = str2num(Array{nEdgeNodeStart+idV});
end

BCele = zeros(numBC,3);
for idBCE = 1:numBC
    BCele(idBCE,:) = str2num(Array{bceStart+idBCE});
end
%meshIEN
%solution_ien(10,:)
end