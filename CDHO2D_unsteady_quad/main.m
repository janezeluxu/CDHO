function [u] = main(fileName,output_filename,solutionFile,vmax)
%This program solve the 2D diffution equation for variable or unifrom order
%simplex elements, the mesh is structured triangles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global itau; 
global Grid_size;
global node_Size;
global TotalDOF;
global a;
global kappa;
%global casenumber;
%global c2;
%global static_factor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[OrderList,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,solution_ien,order_list,edge_usage,...
    iper,uHBCnode,masternodes,slavenodes] ...
    = buildMeshStruct(fileName);
% minA= min(Area);
% maxA = max(Area);
% hmin = minA^0.5;
% alphamin = hmin*norm(a)/(2*norm(kappa));
% 
% hmax = maxA^0.5;
% alphamax = hmax*norm(a)/(2*norm(kappa));
% fprintf('Smallest Mesh Peclet Number: %f\n',alphamin); 
% fprintf('Largest Mesh Peclet Number: %f\n',alphamax); 

%call solveU function, u is the solution field
gGlobal = zeros(TotalDOF,1);
tau= zeros(node_Size,1);
u = zeros(TotalDOF,1);
display('start solveU');

[u,ut] = time_step(vmax,vertexData,BCb,IBC,BCval,meshData,solution_ien,...
    order_list,edge_usage,iper,OrderList,Area,IBC_vertex,uHBCnode);

display('writesolution');
%% write the solution in vtk form
writeMesh(meshData,vertexData,'boundary')
writeSolution(u,solutionFile)
writevtkfile(meshData,vertexData,u(1:node_Size),ut(1:node_Size),tau);
twod_to_vtk_test (output_filename);
end

function [patchArea,patchOrder] = getPatchArea(vertexData,meshData,solution_ien,p_ien,Area)
global node_Size;
patchArea = zeros(node_Size,1);
patchOrder = zeros(node_Size,1);
for node = 1:node_Size 
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
    mineleIEN = zeros(nElement,1);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
        patchArea(node) = patchArea(node)+Area(ele);
        mineleIEN(microEle) = min(pAll);
    end
    patchOrder(node) = min(mineleIEN);
end
end

function [taum] = startingPoint(patchArea,patchOrder,vertexData,IBC)
global node_Size;
global kappa;
global a;
taum = zeros(node_Size,1);
for node = 1:node_Size
    AH = patchArea(node);
    const = 1;  
    p = patchOrder(node);
    PeH = sqrt(const*AH)*norm(a)/(2*norm(kappa)*p);
    %PeHDis(node) = PeH;
    if PeH > 3
        c2 = 1;
        c1 = 1/p/(norm(a));
        %c1=0;
        %c2Distribute(node) = 1;
    else
        %print("diffusion")
        c2 = 2;
        c1 = 1/(p*norm(kappa));
        c1 = 0;
        %c2Distribute(node) = 2;
    end
    
    SurroundEle = vertexData{node,3};
    nEle = length(SurroundEle);
    %const = 1;
    h = (const*AH/nEle)^0.5;
    taum(node) = c1*h^c2;
end
%IBC = IBC_tau{1};
% for i = 1:length(IBC)
%     taum(IBC(i)) = 0;
% end
end