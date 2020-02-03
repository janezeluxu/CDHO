function [Qa] = QaCalculate(OrderList,vertexData,patchArea,patchOrder,...
            IBC_vertex,BCb,meshData,solution_ien,order_list,u,tau,IBC)
global node_Size;
%[~,L1,L2,M1,M2] = solvec1(p,vertexData,faceData,uHBCE,uHBCV,nInt,IENstruct,u);
[taunode,L,M] = solvetauHOmap(OrderList,vertexData,patchArea,patchOrder,...
            IBC_vertex,BCb,meshData,solution_ien,order_list,u);
Qa = zeros(node_Size,1);
for i = 1:node_Size
    Qa(i) = tau(i)*L(i,1)-M(i,1)+tau(i)*L(i,2)-M(i,2);
end

for i = 1:length(IBC)
    j = IBC(i);
    Qa(j) = 0;
end
Qa = (sum(Qa.^2))^0.5;
end