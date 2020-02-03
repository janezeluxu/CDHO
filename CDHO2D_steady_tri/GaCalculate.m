function [Ga,gGlobal] = GaCalculate(p,meshData,solution_ien,...
    order_list,vertexData,BCB,IBC,tau,u,masternodes,slavenodes)
% Assemble of element level matrix   
global TotalDOF; 
global Grid_size;
%global nDirection;
%create uniform p shape function tables, p is a list
[ShapeFuncTable,divSFtable] = ShapeTable(0,p);
gGlobal = zeros(TotalDOF,1);

for ele = 1:Grid_size     
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll),0);
    qPoints = TriGaussPoints(n);
    if range(pAll) == 0
        %disp('uniform p element');
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
    else
        %disp('nonuniform p element');
        %qPoints = TriGaussPoints(n);
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc] = simplexsf.variablePShapeFunctionTable();
        
    end
    
    simplexsfTau = SimplexShapeFunc(qPoints,1);
    [tauShapeFunc, ~] = simplexsfTau.uniformPShapeFunctionTable();
        
    uE = zeros(1,length(IENall));
    for j = 1:length(IENall)
            uE(j) = u(IENall(j));        
    end
    %tauele = tau(ele);
    %qPoints = 0;
    tauele = zeros(3,1);
    for i = 1:3
       tauele(i) = tau(IEN_mesh(i));
    end
    tauELE = max(tauele);
    g = ga(IEN_mesh,IENall,n,vertexData,ShapeFunc,divShapeFunc,tauShapeFunc,tauELE,uE);

    for i = 1:nssl
        rowTemp = IENall(i);
        gGlobal(rowTemp) = g(i) + gGlobal(rowTemp);
    end
      
end

zBc= [IBC,masternodes,slavenodes];
[gGlobal] = strong(zBc,gGlobal);
%[gGlobal] = consistflux(meshData,solution_ien,order_list,vertexData,p,BCB,gGlobal,uE,nDirection);

Ga = (sum(gGlobal.^2))^0.5;
end

function [g] = ga(IEN_mesh,IENall,n,vertexData,ShapeFunc,divShapeFunc,ShapeFunctau,tauele,u)
%element stiffness matrix and force
global kappa;
global a;

quadraturePoints = TriGaussPoints(n);
nP = size(quadraturePoints,1);
sizeN = length(IENall);
g = zeros(sizeN,1);

tri = SimplexGeometry(IEN_mesh,vertexData);
[JInverse, detJ,gij,xCord,yCord,A] =  tri.jacobian();
Jw = detJ*quadraturePoints(:,3);
  
sourceTerm = zeros(nP,1);
for k = 1:nP
    gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
    tau = tauele;
    %SFtau = ShapeFunctau(:,k);
    %tau = SFtau'*tauele;
    g=g+(ShapeFunc(:,k)*(a*(gradNaGlobal'*u')))*Jw(k)...
    +(gradNaGlobal*kappa*(gradNaGlobal'*u'))*Jw(k)...
    +((gradNaGlobal*a')*tau*a*(gradNaGlobal'*u'))*Jw(k)...
     -(ShapeFunc(:,k)*sourceTerm(k))*Jw(k)-(gradNaGlobal*a'*tau*sourceTerm(k))*Jw(k);
    
    
end

end

function [gGlobal] = strong(IBC,gGlobal)
for i = 1:length(IBC)
    j = IBC(i);
    gGlobal(j) = 0;
end
end

function [gGlobal] = consistflux(meshData,solution_ien,order_list,vertexData,p,BCB,gGlobal,u,ndirec)
alpha = 1e-10;
global kappa;
for e = 1:length(BCB)
    bcEle = BCB(e);
    [IEN_mesh,IENall,pAll] = elementIEN(bcEle,meshData,solution_ien,order_list);
    tri = SimplexGeometry(IEN_mesh,vertexData);
    [JInverse, detJ,gij,xCord,yCord,~] =  tri.jacobian();
    %% for each bcb elements
    n = nIntergerPoints(p,0);
    [xi,w] =  GaussQuad(n, 1);
    [edgeNumber,h] = tri.getBCedge(xCord,yCord,alpha);
    qPoints = tri.getBCQuadrature(xi,w,edgeNumber);
    simplexsf = SimplexShapeFunc(qPoints,p);
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
    
    nP = size(qPoints,1);
    sizeN = length(IENall);
    g = zeros(sizeN, 1);    
    Jw = h*qPoints(:,3);
    
    for k = 1:nP
        NaGlobal = ShapeFunc(:,k);
        gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
        g = g - NaGlobal*((kappa*(gradNaGlobal'*u'))'*Jw(k)*ndirec)';
        
    end
    %% Assemble the outlet BC into global matrix
    nssl = length(IENall);
    for i = 1:nssl
        rowTemp = IENall(i);
        gGlobal(rowTemp) = g(i) + gGlobal(rowTemp);
    end
end
end
