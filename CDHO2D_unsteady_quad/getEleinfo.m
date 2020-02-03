function [JInv, detJall,gijall,x,y,A,hA] = getEleinfo(p,meshData,solution_ien,p_ien,vertexData)
global Grid_size;
n = nIntergerPoints(p,0);
quadraturePoints = QuadGaussPoints(n);
nP = size(quadraturePoints,1);
JInv = zeros(Grid_size,nP*4);
detJall = zeros(Grid_size,nP);
gijall = zeros(Grid_size,nP*4);
for ele = 1:Grid_size
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    
    count = 1;
    for k = 1:nP
        
        xi = quadraturePoints(k,1);
        eta = quadraturePoints(k,2);
        
        [JInverse, detJ,gij,x,y,A,hA] =  ...
            getjacobian(IEN_mesh,vertexData,xi,eta);
        
        JInv(ele,count:count+3) = reshape(JInverse,4,1);
        detJall(ele,k) = detJ;
        gijall(ele,count:count+3) = reshape(gij,4,1);
        
        count = count+4;
    end
end
%JIk1 = JI(1,:)
end