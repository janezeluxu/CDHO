function[pOrder,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,solutionIEN,p_list,edge_usage,...
    iper,uHBCnode,masternodes,slavenodes] ...
    = buildMeshStruct(fileName)

global TotalDOF;
global Grid_size;
global node_Size;
global shapeFuncType;
global master;
global slave;
global Dirichlet_edge;
[Grid_size,node_Size,TotalDOF,pOrder,mesh_ien,solution_ien,order_list,...
    edge_usage,edgeNode,vertexDataarray,BC_array,BCb] = readinMesh(fileName);

if strcmp(shapeFuncType, 'H')==1
    edge_usage(edge_usage<=0)=-1;
else
    edge_usage(edge_usage<=0)=1;
end
mesh = geo_Mesh(mesh_ien,vertexDataarray,BC_array,pOrder);
[meshData] = mesh.MeshData();
[vertexData] = mesh.vertexMesh(meshData);
[Area] = mesh.elementArea(meshData,vertexData);
[IBC,BCval] = mesh.BoundaryCondition();
%IBC = BC_array;
%BCval = zeros(length(IBC),1);
ien = IEN(mesh_ien,solution_ien,order_list,0);
[solutionIEN,p_list] = ien.Construct_IEN();

IBC_vertex = [];
for i = 1:length(IBC)
    if IBC(i)<node_Size
        IBC_vertex = [IBC_vertex,IBC(i)];
    end
end
%uHBCnode_u = findtaunode(BCb,mesh_ien,IBC,vertexData,BCval);
%uHBCnode_v = findtaunode(BCb,mesh_ien,IBC,vertexData,BCval);
uHBCnode = findtaunode(BCb,mesh_ien,IBC,vertexData,BCval,Dirichlet_edge);
%uHBCnode = [];
%Build periodic
iper = zeros(TotalDOF,1); 
masternodes = [];
slavenodes = [];
% [masternodes,slavenodes]=getiper(edgeNode,master,slave,BC_array);
% for i = 1:TotalDOF
%     iper(i) = i;
%     for j = 1:length(masternodes)
%         iper(slavenodes(j)) = masternodes(j);
%     end       
% end
end

function uHBCnode = findtaunode(BCb,ien_map,IBC,vertexData,BCval,Dirichlet_edge)
global node_Size;
uHBCnode = zeros(node_Size,3);
nodeList = [];
%% get all boundary element node
for i = 1:length(BCb)
    ele = BCb(i,1);
    edgeGeoClassify = BCb(i,2);
    if ismember(edgeGeoClassify,Dirichlet_edge) 
        IEN_mesh = ien_map(ele,:);
        nodeList = [nodeList,IEN_mesh];
    end
end
nodeList = setdiff(nodeList,IBC);
%% set uHBCnode flag
for i = 1:length(nodeList)
    bcnode = nodeList(i);
    uHBCnode(bcnode,1) = 1;
    [fitnode] = fitExactV(vertexData,bcnode,IBC);
    uHBCnode(bcnode,2) = fitnode;
    [Lia,loc] = ismember(fitnode,IBC);
    uHBCnode(bcnode,3) = BCval(loc);
end
%uHBCnode %= zeros(node_Size,3);
end

