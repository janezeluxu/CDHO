function [concentration]=LinePlot(fileName,solution,x,y)
[OrderList,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,uHBCE,solution_ien,order_list,edge_usage,iper] ...
    = buildMeshStruct(fileName);
 u = load(solution);
 
 concentration = zeros(1,length(y));
 
 %% start the search from boundary elements
 [startingEle]=findNextEle(x(1),y(1),BCb,...
        meshData,solution_ien,order_list,vertexData)
 %startingEle = 35;
 lastEle = startingEle;
 
 for i = 1:length(x)
    eleSearch = getSurrondEle(lastEle,meshData,...
        solution_ien,order_list,vertexData);
    [newEle] = findNextEle(x(i),y(i),eleSearch,...
        meshData,solution_ien,order_list,vertexData);
    
    concentration(i) = getSolution(u,x(i),y(i),newEle,meshData,...
        solution_ien,order_list,vertexData);
    lastEle = newEle;
 end
 lastEle
 %plot(y,concentration,'b*-')
end

function concentration = getSolution(u,x,y,Ele,meshData,solution_ien,p_ien,vertexData)
[ien_mesh,IENall,pAll] = elementIEN(Ele,meshData,solution_ien,p_ien);
tri = SimplexGeometry(ien_mesh,vertexData);
[JInverse, detJ,gij,xCord,yCord,A,hA] =  tri.jacobian();
EvalPoints = getPoints(x,y,xCord,yCord);
simplexsf = SimplexShapeFunc(EvalPoints,pAll(1));
[ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable(); 

solution = zeros(length(IENall),1);
for i = 1:length(IENall)
    solution(i) = u(IENall(i));
end
concentration = solution'*ShapeFunc;
end

function eleSearch = getSurrondEle(lastEle,meshData,solution_ien,p_ien,vertexData)
[ien_mesh] = elementIEN(lastEle,meshData,solution_ien,p_ien);
eleSearch = [];
for i = 1:3
    SurroundEle = vertexData{ien_mesh(i),3};
    eleSearch = [eleSearch,SurroundEle];
end
end

function [searchELE]=findNextEle(x,y,searchPatch,meshData,solution_ien,p_ien,vertexData)
%% search for the element contains x y point
searchELE = 0;
for i = 1:length(searchPatch)
    ele = searchPatch(i);
    [ien_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    tri = SimplexGeometry(ien_mesh,vertexData);
    [JInverse, detJ,gij,xCord,yCord,A,hA] =  tri.jacobian();
    [evaluatePoints,A1,A2,A3]=getPoints(x,y,xCord,yCord);
    lam2 = evaluatePoints(1);
    lam3 = evaluatePoints(2);
    lam1 = 1-lam2-lam3;
    if (lam2-1<=1e-10 && lam2>=-1e-10 && lam3>=-1e-10 && lam3-1<=1e-10...
            && lam1>=-1e-10 && lam1-1<=1e-10)
        searchELE = ele;
        
    end
end
end


function [evaluatePoints,A1,A2,A3]=getPoints(x,y,xCordM1,yCordM1)
%map to area quadinates
evaluatePoints = zeros(1,2);
%for i = 1:length(x)
i = 1;
    x1 = xCordM1(1);
    x2 = xCordM1(2);
    x3 = xCordM1(3);
    y1 = yCordM1(1);
    y2 = yCordM1(2);
    y3 = yCordM1(3);
    Atotal = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = x(i);
    x2 = xCordM1(2);
    x3 = xCordM1(3);
    y1 = y(i);
    y2 = yCordM1(2);
    y3 = yCordM1(3);
    A1 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = xCordM1(1);
    x2 = x(i);
    x3 = xCordM1(3);
    y1 = yCordM1(1);
    y2 = y(i);
    y3 = yCordM1(3);
    A2 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = xCordM1(1);
    x2 = xCordM1(2);
    x3 = x(i);
    y1 = yCordM1(1);
    y2 = yCordM1(2);
    y3 = y(i);
    A3 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    lam2 = A2/Atotal;
    lam3 = A3/Atotal;
    evaluatePoints(i,:) = [lam2,lam3];
%end
%% check if A1/Atotal = 1-lam2-lam3
A1test = A1/Atotal;
lam1 = 1-lam2-lam3;
end