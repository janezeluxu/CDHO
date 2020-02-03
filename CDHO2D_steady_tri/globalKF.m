function [kGlobal, fGlobal] = globalKF(p,meshData,solution_ien,p_ien,edge_usage,BCb,IBC,BCval,...
                                vertexData,iper,tau)
% Assemble of element level matrix   
global TotalDOF; 
global Grid_size;
%create uniform p shape function tables, p is a list
[ShapeFuncTable,divSFtable] = ShapeTable(0,p);

row = [];
column = [];
value = [];
fGlobal = zeros(TotalDOF,1);
%Kindex = ones(TotalDOF,TotalDOF)*-1;
%Kindex = spalloc(TotalDOF,TotalDOF,TotalDOF*200);
keyCount = 0;

for ele = 1:Grid_size 
%for ele = 1:2    
    %ele
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions    
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    
    edgeU = edge_usage(ele,:);
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll),0);
    if range(pAll) == 0
        %disp('uniform p element');
        ShapeFunc = ShapeFuncTable{pAll(1)};
        divShapeFunc = divSFtable{pAll(1)};
        
        %simplexsf = SimplexShapeFunc(qPoints,pAll(1));
        %[ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();  
    else
        disp('nonuniform p element');
        qPoints = TriGaussPoints(n);
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc] = simplexsf.variablePShapeFunctionTable();   
    end

%     count = 3;
%     for i = 1:3
%         for j = 1:p-1
%             count = count+1;
%             ShapeFunc(count,:)= edgeU(i)*ShapeFunc(count,:);
%             divShapeFunc(count,:,:)= edgeU(i)*divShapeFunc(count,:,:);
%         end
%     end
            
    %find element Shape Functions, insert stiffness matrix and force vectors
    %into global matrixs(spars)
    %tauele = tau(ele);
    
    tauele = zeros(3,1);
    for i = 1:3
       tauele(i) = tau(IEN_mesh(i));
    end
    tauELE = max(tauele);
    qPoints = 0;
    
    simplex = SimplexStiffMatrix(ele,IEN_mesh,IENall,pAll,n,qPoints,vertexData,...
        ShapeFunc,divShapeFunc,tauELE);
    %ShapeFunc
    %divShapeFunc
    [elementK,elementF] = simplex.eleStiffMatrix();
    
    %sparse matrix input row, column, and value list, using row and column
    %combination as key
    for i = 1:nssl
        rowTemp = IENall(i);
        fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
        for j = 1:nssl           
            columnTemp = IENall(j);
            %from a given key, find the value index, -1 means it is a new
            %key, thus insert rowTemp, columnTemp, value at the end of the
            %list
            %indexKglobal = Kindex(rowTemp,columnTemp);
            %indexKglobal = IENindex(rowTemp, columnTemp, row, column);
            %if indexKglobal == 0
            %if indexKglobal == -1
                keyCount = keyCount+1; %its a new key
                row = [row, rowTemp]; 
                column = [column, columnTemp]; 
                value = [value,elementK(i,j)]; 
                %Kindex(rowTemp,columnTemp) = keyCount; %update Kindex
            %else
                %value(indexKglobal) = elementK(i,j)+value(indexKglobal);
            %end
        end
    end
    
%     for i = 1:nssl
%         rowTemp = IENall(i);
%         fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
%         for j = 1:nssl           
%             columnTemp = IENall(j);
%             %from a given key, find the value index, -1 means it is a new
%             %key, thus insert rowTemp, columnTemp, value at the end of the
%             %list
%             indexKglobal = Kindex(rowTemp,columnTemp);
%             %indexKglobal = IENindex(rowTemp, columnTemp, row, column);
%             if indexKglobal == 0
%             %if indexKglobal == -1
%                 keyCount = keyCount+1; %its a new key
%                 row = [row, rowTemp]; 
%                 column = [column, columnTemp]; 
%                 value = [value,elementK(i,j)]; 
%                 Kindex(rowTemp,columnTemp) = keyCount; %update Kindex
%             else
%                 value(indexKglobal) = elementK(i,j)+value(indexKglobal);
%             end
%         end
%     end
      
end

display('nowBC');
boundary = BC(BCb,IBC,BCval,iper,row,column,value,fGlobal,...
    vertexData,p,solution_ien,meshData,p_ien,tau);
[row,column,value,fGlobal] = boundary.ApplyBC();

% [row,column,value,fGlobal] = StrongBC(IBC,BCval,...
%                 row,column,value,fGlobal);
 kGlobal = sparse(row, column, value);
end

