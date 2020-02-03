function [tau,c1] = solvec1(p,vertexData,faceData,IBC,uHBCE,nInt,meshData,solution_ien,order_list,u)
global c2;
global kappa; 
global a;
global node_Size;
global Grid_size;
[ShapeFuncTable,divSFtable] = ShapeTable(nInt,p);
[SFTable,divSF] = ShapeTable(nInt,p);
c1 = zeros(node_Size,1);

for node = 1:node_Size
    LS = LeastSquare(node,vertexData,u,uHBCE);
    
    a0 = LS(1);
    a1 = LS(2);
    a2 = LS(3);
      
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
    Area = zeros(1,nElement);
    xyc = zeros(nElement,3);
    xyc(:,1) = 1;
       
    AA = [0;0];
    EE = [0;0];
    FF = [0;0];
    for microEle = 1:nElement
        %[IENall,pAll] = elementIEN(SurroundEle(microEle),IENstruct,p);
        ele = SurroundEle(microEle);
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
        %c1Ele = 0;
        %simplex = SimplexStiffMatrix(SurroundEle(microEle),IENstruct,IENall,nInt,0,vertexData,ShapeFuncTable,divSFtable,c1Ele);
        %[JInverse, detJ,~,xCord,yCord] =  simplex.jacobian();
        Area(microEle) = faceData(SurroundEle(microEle));
        
        tri = SimplexGeometry(IEN_mesh,vertexData);
        
        [JInverse, detJ,gij,xCord,yCord,A] =  tri.jacobian();
        he = tri.hCalculate(gij,Area(microEle));
        %Area(microEle)=A;
        %Aeror = faceData{SurroundEle(microEle),4}-A
        %----------------------------- uh ----------------%
        
        uA = zeros(1,length(IENall));
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
        end
        uh = uA*ShapeFunc;
        duh = uA*divShapeFunc*JInverse;
        
        %he = (Area(microEle))^0.5;
        
        AA = AA+(a'*a)*he^(c2)*duh'*detJ;
        EE = EE+a'*uh*detJ;
        FF = FF + kappa*duh'*detJ;
        

        SF = SFTable{1};
        x = SF'*xCord';
        y = SF'*yCord';
        xyc(microEle,2) = x;
        xyc(microEle,3) = y;
    end
    uc = xyc*LS;
    uH = Area*uc;    
    sumA = sum(Area);    

    He = sumA^0.5;
    BB = He^(c2)*(a'*a)*[a1;a2]*sumA;
    CC = -a'*uH;
    DD = kappa*[a1;a2]*sumA;
    

    M = AA-BB;
    L = CC+DD+EE-FF;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c1(node) = M\L;
     
     if c1(node) < 0
         c1(node)=0;
         %c1(node)=abs(c1(node));
     end
end
tau= zeros(Grid_size,1);
for ele = 1:Grid_size
    Area = faceData(ele);
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
    tri = SimplexGeometry(IEN_mesh,vertexData);
    %tri = SimplexGeometry(ele,IENstruct,vertexData);
    [JInverse, detJ,gij,xCord,yCord,A] =  tri.jacobian();
    he = tri.hCalculate(gij,Area);
    tau1 = c1(IENall(1))*he^c2;
    tau2 = c1(IENall(2))*he^c2;
    tau3 = c1(IENall(3))*he^c2;
    tau(ele) = max([tau1,tau2,tau3]);
end
for i = 1:length(IBC)
    j = IBC(i);
    tau(j) = 0;
    c1(j) = 0;
end
end

function LS = LeastSquare(node,vertexData,u,uHBCE)
    global BCType;
    SurroundNode = vertexData{node,4}; %sorrounding node
    nNode = length(SurroundNode);
    CoordMatrix = zeros(nNode,3);
    uE = zeros(nNode,1);
    CoordMatrix(:,1) = 1;
    for micronode = 1:nNode        
        CoordMatrix(micronode,2:3)= vertexData{SurroundNode(micronode),2};
        uE(micronode) = u(SurroundNode(micronode));
    end
    LScoef = CoordMatrix\uE;
    %% strongly impose bc
    %if (strcmp(BCType{1},'Strong')==1||strcmp(BCType{1},'flux')==1)
        %disp ('Strong!!!');
        if any(node ==uHBCE(1,:)) ==1
            a1 = LScoef(2);
            a2 = LScoef(3);
            
            index = find(node ==uHBCE(1,:));
            Cord = vertexData{uHBCE(2,index),2};
            index = index(1);
            Value = u(uHBCE(2,index));
            a0 = Value-a1*Cord(1)-a2*Cord(2);
            LS  = [a0;a1;a2];
        else
            LS = LScoef;
        end
    %end
end