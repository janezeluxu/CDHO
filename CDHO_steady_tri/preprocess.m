function[OrderList,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,uHBCE,solution_ien,order_ien,edge_usage,iper] ...
    = preprocess(fileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global kappa; 
global a;
global Grid_size;
global casenumber;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the mesh data structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[OrderList,meshData,vertexData,Area,IBC,IBC_vertex,BCval,BCb,uHBCE,...
    solution_ien,order_ien,edge_usage,iper] ...
    = buildMeshStruct(fileName);
minA= min(Area);
maxA = max(Area);
hmin = minA^0.5;
alphamin = hmin*norm(a)/(2*norm(kappa));

hmax = maxA^0.5;
alphamax = hmax*norm(a)/(2*norm(kappa));
OrderList
fprintf('Smallest Mesh Peclet Number: %f\n',alphamin); 
fprintf('Largest Mesh Peclet Number: %f\n',alphamax); 

if casenumber ==100
    hA = zeros(Grid_size,1);
    for e = 1:Grid_size
        [IEN_mesh,IENall,pAll] = elementIEN(e,meshData,solution_ien,order_ien);
        tri = SimplexGeometry(IEN_mesh,vertexData);
        [JInverse, detJ,gij,xCord,yCord,A,hA(e)] =  tri.jacobian();
    end 
    mean(hA)
end