function [kGlobal, fGlobal] = globalKF(p,meshData,solution_ien,p_ien,edge_usage,BCb,IBC,BCval,...
                                vertexData,iper,tau)
% Assemble of element level matrix   
global TotalDOF; 
global Grid_size;
%create uniform p shape function tables, p is a list
n = nIntergerPoints(p,0);
qPoints = QuadGaussPoints(n);

[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);

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
    
    nssl = length(IENall);
    %n = nIntergerPoints(max(pAll),0);
    
    %find element Shape Functions, insert stiffness matrix and force vectors
    %into global matrixs(spars)
    %tauele = tau(ele);
    
    tauele = zeros(3,1);
    for i = 1:3
       tauele(i) = tau(IEN_mesh(i));
    end
    tauELE = max(tauele);
    %qPoints = 0;
    quad = quadStiffMatrix(ele,IEN_mesh,IENall,pAll,n,qPoints,vertexData,...
        ShapeFunc,divShapeFunc,tauELE);
    [elementK,elementF] = quad.eleStiffMatrix();
    
    %sparse matrix input row, column, and value list, using row and column
    %combination as key
    for i = 1:nssl
        rowTemp = IENall(i);
        fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
        for j = 1:nssl
            columnTemp = IENall(j);
            keyCount = keyCount+1; %its a new key
            row = [row, rowTemp];
            column = [column, columnTemp];
            value = [value,elementK(i,j)];
        end
    end
      
 end

display('nowBC');
boundary = BC(BCb,IBC,BCval,iper,row,column,value,fGlobal,...
    vertexData,p,solution_ien,meshData,p_ien,tau);
[row,column,value,fGlobal] = boundary.ApplyBC();

% [row,column,value,fGlobal] = StrongBC(IBC,BCval,...
%                 row,column,value,fGlobal);
 kGlobal = sparse(row, column, value);
end


function [row,column,value,fGlobal] = StrongBC(IBC,BCval,...
                row,column,value,fGlobal)
            
            % Account for BC in forcing vector.
            
            for i = 1:length(column)
                [BCvalIndex,a,BCvalIndex]=intersect(column(i),IBC);
                if BCvalIndex~=0
                    fGlobal(row(i))= fGlobal(row(i))-value(i)*BCval(BCvalIndex);
                end
            end
            
            %add Boundary conditions from IBC
            for i = 1:length(row)
                if (any(row(i) == IBC))
                    row(i) = -1;
                    column(i) = -1;
                    value(i) = NaN;
                end
                if (any(column(i) == IBC))
                    row(i) = -1;
                    column(i) = -1;
                    value(i) = NaN;
                end
            end
            
            %add back value as 1 at [IBC(i),IBC(i)]
            for i = 1: length (IBC)
                rowTemp = IBC(i);
                fGlobal(rowTemp) = BCval(i);
                %fGlobal(rowTemp) = 0;
                row = [row,rowTemp];
                column = [column,rowTemp];
                value = [value,1];
            end
            
            
            row = row(row~=-1);
            column = column(column~=-1);
            value = value(~isnan(value));
            
        end

