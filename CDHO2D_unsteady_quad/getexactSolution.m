function [u] = getexactSolution(maxIter,fileName,output_filename,solutionFile,GaTol)
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
    IBC,IBC_vertex,BCval,BCb,uHBCnode,solution_ien,order_list,edge_usage,...
    iper,masternodes,slavenodes] ...
    = buildMeshStruct(fileName);
minA= min(Area);
maxA = max(Area);
hmin = minA^0.5;
alphamin = hmin*norm(a)/(2*norm(kappa));

hmax = maxA^0.5;
alphamax = hmax*norm(a)/(2*norm(kappa));
fprintf('Smallest Mesh Peclet Number: %f\n',alphamin); 
fprintf('Largest Mesh Peclet Number: %f\n',alphamax); 

%call solveU function, u is the solution field
gGlobal = zeros(TotalDOF,1);
tau= zeros(node_Size,1);
u = zeros(TotalDOF,1);
display('start solveU');

[u] = solveU(vertexData,BCb,IBC,BCval,meshData,solution_ien,order_list,...
     edge_usage,iper,tau,OrderList);
end