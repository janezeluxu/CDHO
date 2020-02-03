function [masternodes,slavenodes]=getiper(edgeNode,master,slave,BC_array)

%edgeNode = getedgeNode(fileName);

master_n = findnode(master(1),BC_array);
slave_n = findnode(slave(1),BC_array);
master_nm1 = master_n;
slave_nm1 = slave_n;

master_other = findnode(master(3),BC_array);
slave_other = findnode(slave(3),BC_array);

[masterstop]=findedge(master_other,master_other,edgeNode,master(2),BC_array);

masternodes = [];
slavenodes = [];

for i = 1:100
    
[master_np1,allnodemaster,usagemaster,slave_np1,allnodeslave,usageslave]...
    = periodicnode(master_n,slave_n,master_nm1,slave_nm1,BC_array,edgeNode,master,slave);

[masternodes,slavenodes]=buildperiodic(allnodemaster,allnodeslave,...
    usagemaster,usageslave,masternodes,slavenodes);

master_nm1 = master_n;
slave_nm1 = slave_n;
master_n = master_np1;
slave_n = slave_np1;

if master_np1 == masterstop
    break
end
end
masternodes = [masternodes,master_other];
slavenodes = [slavenodes,slave_other];
end

function [opponodemaster,allnodemaster,usagemaster,opponodeslave,allnodeslave,usageslave]...
    = periodicnode(lastmaster,lastslave,previousmaster,previousslave,BC_array,edgeNode,master,slave)

%startingmaster = findnode(master(1),BC_array);
[opponodemaster,allnodemaster,usagemaster] = findedge(lastmaster,previousmaster,edgeNode,master(2),BC_array);
%startingslave = findnode(slave(1),BC_array);
[opponodeslave,allnodeslave,usageslave] = findedge(lastslave,previousslave,edgeNode,slave(2),BC_array);

end

function [masternodes,slavenodes]=buildperiodic(allnodemaster,allnodeslave,...
    usagemaster,usageslave,masternodes,slavenodes)
%masternodes = [];
if usagemaster ==0
    edgenodemaster = [allnodemaster(2),allnodemaster(1),flip(allnodemaster(3:end))];
    masternodes = [masternodes,edgenodemaster(2:end)];
elseif usagemaster ==1
    masternodes = [masternodes,allnodemaster(2:end)];
end
%slavenodes = [];
if usageslave ==0
    edgenodeslave = [allnodeslave(2),allnodeslave(1),flip(allnodeslave(3:end))];
    slavenodes = [slavenodes,edgenodeslave(2:end)];
elseif usageslave ==1
    slavenodes = [slavenodes,allnodeslave(2:end)];
end

end

function [opponode,allnode,usage]=findedge(node,previousNode,edgeNode,geomClassi,BC_array)
%% find node belonged edge
index1 = find(node==edgeNode(:,1));
index2 = find(node==edgeNode(:,2));
count = 1;
for i = 1:length(index1)
    edge(count,:) = edgeNode(index1(i),:);
    count = count+1;
end
for i = 1:length(index2)
    edge(count,:) = edgeNode(index2(i),:);
    count = count+1;
end


for i = 1:size(edge,1)
    if edge(i,1)==node
        if BC_array(edge(i,2))==geomClassi
            if edge(i,2) ~=previousNode
                opponode = edge(i,2);
                allnode = edge(i,:);
                usage = 1;
            end
        end
    elseif edge(i,2)==node
        if BC_array(edge(i,1))==geomClassi
            if edge(i,1) ~=previousNode
                allnode = edge(i,:);
                opponode = edge(i,1);
                usage = 0;
            end
        end
    end
end

end

function node = findnode(goem,IDall)
node = [];
for i = 1:size(IDall,1)
    nodeID = i;
    node_tag = IDall(i);
    for j = 1:size(goem,1)
        if (node_tag == goem)
            node = [node,nodeID];
        end
    end
end
end

function edgeNode = getedgeNode(fileName)
fileID = fopen(fileName,'r');
tline = fgetl(fileID);
i = 1;

while ischar(tline)
    C = textscan(tline,'%s');
    c= C{1};
    
    if ( strcmp(c{1},'edge_Size')==1)
        edgeNum = str2num(c{2});
    elseif ( strcmp(c{1},'Solution_Order')==1)
        npstart = i;
    end
    
    if ( strcmp(c{1},'EdgeNode')==1)
        nEdgeNodeStart = i;
    end

    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);
maxOrder = str2num(Array{npstart+1});
%nodeSizet = 83;
edgeNode = zeros(edgeNum,maxOrder+1);
for idV = 1:edgeNum
    edgeNode(idV,:) = str2num(Array{nEdgeNodeStart+idV});
end
end