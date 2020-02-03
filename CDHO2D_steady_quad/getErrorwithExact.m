function getErrorwithExact(fileName,solution)
global TotalNode;
[OrderList,meshData,vertexData,Area,IBC,IBC_V,BCval,BCb,uHBCE,...
    solution_ien,order_ien,edge_usage,iper] ...
    = buildMeshStruct(fileName);

 variable=load(solution);
 nInt = 200;
 [L2,H1] = errorEstimate(variable,vertexData,meshData,solution_ien,order_ien,...
     nInt,OrderList)
end