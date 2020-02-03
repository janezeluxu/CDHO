function [concentration]=LinePlot(fileName,solution,x,y)

[OrderList,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,uHBCnode,solution_ien,order_list,edge_usage,...
    iper,masternodes,slavenodes] ...
    = buildMeshStruct(fileName);

 u = load(solution);
 
 concentration = zeros(1,length(y));
 
 %% start the search from boundary elements
 [startingEle]=findNextEle(x(1),y(1),BCb,...
        meshData,solution_ien,order_list,vertexData);
 startingEle = 2;
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
 
 %plot(y,concentration,'b*-')
end

function concentration = getSolution(u,x,y,Ele,meshData,solution_ien,p_ien,vertexData)
[ien_mesh,IENall,pAll] = elementIEN(Ele,meshData,solution_ien,p_ien);
[JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                    getjacobian(ien_mesh,vertexData,0,0);
EvalPoints = getPoints(x,y,xCord,yCord);

[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(EvalPoints,pAll(1));

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
distanceToCenter = zeros(length(searchPatch),1);
for i = 1:length(searchPatch)
    ele = searchPatch(i);
    [ien_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);

    [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                    getjacobian(ien_mesh,vertexData,0,0);
                
    [EvalPoints, distanceToCenter(i)]=getPoints(x,y,xCord,yCord);
end
[value,index]=min(distanceToCenter);
searchELE = searchPatch(index);
end


function [EvalPoints,distanceToCenter]=getPoints(x,y,xCordM1,yCordM1)
%% get center coordinate 
x1 = xCordM1(1);
x2 = xCordM1(2);
x3 = xCordM1(3);
x4 = xCordM1(4);

y1 = yCordM1(1);
y2 = yCordM1(2);
y3 = yCordM1(3);
y4 = yCordM1(4);

xmid = (x1+x2+x3+x4)/4;
ymid = (y1+y2+y3+y4)/4;

distanceToCenter = sqrt((x-xmid)^2+(y-ymid)^2);

xi = abs(x-x1)/abs(x1-x2);
eta = abs(y-y1)/abs(y4-y1);
EvalPoints(1) = xi;
EvalPoints(2) = eta;

% %% get evalPoints in xi,eta coordinate
% a = (y4-y1)/(x4-x1);
% b = y1-a*x1;
% l1 = abs(a*x-y+b)/sqrt(a^2+1);
% l = sqrt((x1-x2)^2+(y2-y1)^2);
% xi = l1/l;
% a = (y2-y1)/(x2-x1);
% b = y1-a*x1;
% m1 = abs(a*x-y+b)/sqrt(a^2+1);
% m = sqrt((x1-x4)^2+(y4-y1)^2);
% eta = m1/m;
% EvalPoints(1) = xi;
% EvalPoints(2) = eta;
end