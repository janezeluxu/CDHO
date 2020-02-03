classdef BC < handle
    
    properties (Access = private)
        bcb;
        ibc;
        bcval;
        IPER;
        row;
        column;
        value;
        fglobal;
        edge;
        P;
        ien;
        ien_mesh;
        order_ien;
        tau;
        VERTEX;
    end
    
    methods(Access = public)
        
        function BoundaryCondition = BC(bcb,ibc,bcval,Iper,Row,Column,Value,...
                FGlobal,vertex,Order,solution_ien,mesh_map,p_ien,Tau)
            BoundaryCondition.bcb = bcb;
            BoundaryCondition.ibc = ibc;
            BoundaryCondition.bcval = bcval;
            BoundaryCondition.IPER = Iper;
            BoundaryCondition.row = Row;
            BoundaryCondition.column = Column;
            BoundaryCondition.value = Value;
            BoundaryCondition.fglobal = FGlobal;
            
            BoundaryCondition.VERTEX = vertex;
            BoundaryCondition.P = Order;
            BoundaryCondition.ien = solution_ien;
            BoundaryCondition.ien_mesh = mesh_map;
            BoundaryCondition.order_ien=p_ien;
            BoundaryCondition.tau = Tau;
        end
        
        function [row,column,value,fGlobal] = ApplyBC(BoundaryCondition)
            %global nDirection;
            row = BoundaryCondition.row;
            column = BoundaryCondition.column;
            value = BoundaryCondition.value;
            fGlobal = BoundaryCondition.fglobal;
            
            BCB = BoundaryCondition.bcb;
            IBC = BoundaryCondition.ibc;
            BCval = BoundaryCondition.bcval;
            
            [row,column,value,fGlobal] = BoundaryCondition.StrongBC(IBC,BCval,...
                        row,column,value,fGlobal);
            %[row,column,value,fGlobal] = BoundaryCondition.consistflux(BCB,...
             %           row,column,value,fGlobal,nDirection);
             %[row,column,value,fGlobal] = BoundaryCondition.PeriodicBC(row,column,...
             %    value,fGlobal);
            
        end
        
        
        function [row,column,value,fGlobal] = StrongBC(BoundaryCondition,IBC,BCval,...
                row,column,value,fGlobal)
            
            % Account for BC in forcing vector.
            global TotalDOF
            
%             for i = 1:length(IBC)
%                 for j = 1:TotalDOF
%                     % XXX use ien instead of ismember and TotalDOF, will be faster
%                     if(~ismember(j,IBC) || BCval(i)==0)    % Don't do this if there's a BC here, it'll be 0 anyways
%                         tmpIndex = BoundaryCondition.IENindex(j,IBC(i),row,column);
%                         %tmpIndex = Kindex(j,IBC(i));
%                         if(tmpIndex>0)
%                             fGlobal(j) = fGlobal(j) - BCval(i)*value(tmpIndex);
%                         end
%                     end
%                 end
%             end
            
            BCIndex = zeros(TotalDOF,1);
            for j = 1:TotalDOF
                [IBCval,a,BCvalI]=intersect(j,IBC);
                if BCvalI~=0
                    BCIndex(j)=BCvalI;
                end
            end
            
            for i = 1:length(column)
                %[IBCvalue,a,BCvalIndex]=intersect(column(i),IBC);
                c = column(i);
                IC = BCIndex(c);
                if IC~=0
                    fGlobal(row(i))= fGlobal(row(i))-value(i)*BCval(IC);
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
              
        function [row,column,value,fGlobal] = PeriodicBC(BoundaryCondition,row,column,value,fGlobal)
            iper = BoundaryCondition.IPER;
            global TotalDOF
            iperBC = [];
            % Account for periodicity
            for i = 1:length(iper)
                if(iper(i) ~= i)
                    rowTemp = iper(i);
                    iperBC = [iperBC,i];
                    %rowTemp is master, i is slave
                    % forcing vector
                    fGlobal(rowTemp) = fGlobal(rowTemp) + fGlobal(i);
                    fGlobal(i) = 0;
                    
                    % add slave row to master row in stiffness matrix
                    for j = 1:TotalDOF
                        %IndexSlave = Kindex(i,j);
                        IndexSlave = BoundaryCondition.IENindex(i,j,row,column);
                        %IndexMaster = Kindex(rowTemp,j);
                        %IndexMaster = BoundaryCondition.IENindex(rowTemp,j,row,column);
                        columnPrim = iper(j);
                        
                        if( IndexSlave~=0)
                            for ll = 1:length(IndexSlave)
                                row = [row, rowTemp];
                                column = [column, columnPrim];
                                value = [value,value(IndexSlave(ll))];
                            end
                        end
                        
                    end
                end
            end
            
            %% remove slave nodes from iperBC
            for i = 1:length(row)
                if  (any(row(i) == iperBC))
                    row(i) = -1;
                    column(i) = -1;
                    value(i) = NaN;
                end
            end
            %fGlobal
            
            %% add slave == master node
            for i =1:length(iperBC)
                row = [row,iperBC(i)];
                column = [column,iperBC(i)];
                value = [value,1];
                
                row = [row,iperBC(i)];
                column = [column,iper(iperBC(i))];
                value = [value,-1];
            end
            
            row = row(row~=-1);
            column = column(column~=-1);
            value = value(~isnan(value));
            
        end
            
        function [row,column,value,fGlobal] = flux(BoundaryCondition,row,column,value,fGlobal)
            
            row = row;
            column = column;
            value = value;
            fGlobal = fGlobal;
        end
        
        function [row,column,value,fGlobal] = consistflux(BoundaryCondition,bcb,row,column,value,fGlobal,ndirec)
            solution_ien = BoundaryCondition.ien;
            meshData = BoundaryCondition.ien_mesh;
            vertexData = BoundaryCondition.VERTEX;
            p = BoundaryCondition.P;
            
            alpha = 1e-10;
            for e = 1:length(bcb)
                bcEle = bcb(e);
                [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
                tri = SimplexGeometry(IEN_mesh,vertexData);
                [JInverse, detJ,gij,xCord,yCord,~] =  tri.jacobian();
                %% for each bcb elements
                n = nIntergerPoints(p,0);
                [xi,w] =  GaussQuad(n, 1);
                [edgeNumber,h] = tri.getBCedge(xCord,yCord,alpha);
                qPoints = tri.getBCQuadrature(xi,w,edgeNumber);
                simplexsf = SimplexShapeFunc(qPoints,p);
                [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
                
                simplex = SimplexStiffMatrix(bcEle,IEN_mesh,IENall,n,qPoints,vertexData,...
                    ShapeFunc,divShapeFunc,0);
                [elementK,~] = simplex.fluxBcStiffMatrix(ndirec,h);
                
                nssl = length(IENall);
                %% Assemble the outlet BC into global matrix
                for i = 1:nssl
                    rowTemp = IENall(i);
                    for j = 1:nssl
                        columnTemp = IENall(j);
                        indexKglobal = BC.IENindex(rowTemp, columnTemp, row, column);
                        if indexKglobal == 0
                            row = [row, rowTemp];
                            column = [column, columnTemp];
                            value = [value,elementK(i,j)];
                        else
                            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
                        end
                    end
                end
            end
            fGlobal = fGlobal;
        end

        function [row,column,value,fGlobal] = Inlet(BoundaryCondition,BCval,row,column,value,fGlobal,InletBC,InletEdge,nin)
            IENstruct = BoundaryCondition.ienstruct;
            p = BoundaryCondition.P;
            c1 = BoundaryCondition.C1;
            vertexData = BoundaryCondition.VERTEX;
            edgeD = BoundaryCondition.edge;
            %% calculate added Stiffness matrix for inlet
            for ii = 1:length(InletBC)
                ele = InletBC(ii);
                BCvalEle = BCval(ii);
                [IENall,pAll] = elementIEN(ele,IENstruct);
                vIDs = edgeD{InletEdge(ii),2};
                [~,edgeNumber] = find(vIDs(1) == IENall);
                
                n = nIntergerPoints(p,0);
                [xi, w] = GaussQuad(n, 1);
                qPoints = BC.getBCQuadrature(xi,w,edgeNumber);
                
                simplexsf = SimplexShapeFunc(qPoints,p);
                [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
                
                nssl = length(IENall);
                c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
                simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,qPoints,vertexData,ShapeFunc,divShapeFunc,c1Ele);
                
                [elementK,elementF] = simplex.WeakBCInletStiffMatrix(nin,BCvalEle);
                %Assemble the inlet BC into global matrix
                for i = 1:nssl
                    rowTemp = IENall(i);
                    fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
                    for j = 1:nssl
                        columnTemp = IENall(j);
                        indexKglobal = BC.IENindex(rowTemp, columnTemp, row, column);
                        if indexKglobal == 0
                            %if indexKglobal == -1
                            %keyCount = keyCount+1; %its a new key
                            row = [row, rowTemp];
                            column = [column, columnTemp];
                            value = [value,elementK(i,j)];
                            %Kindex(rowTemp,columnTemp) = keyCount; %update Kindex
                        else
                            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
                        end
                    end
                end
            end
            
        end
               
        function [row,column,value,fGlobal] = Outlet(BoundaryCondition,BCval,row,column,value,fGlobal,InletBC,OutEdge,nout)
            IENstruct = BoundaryCondition.ienstruct;
            p = BoundaryCondition.P;
            c1 = BoundaryCondition.C1;
            vertexData = BoundaryCondition.VERTEX;
            edgeD = BoundaryCondition.edge;
            
            %% calculate added Stiffness matrix for outlet
            for ii = 1:length(InletBC)
                ele = InletBC(ii);
                BCvalEle = BCval(ii);
                %BCvalEle = 0;
                [IENall,pAll] = elementIEN(ele,IENstruct);
                vIDs = edgeD{OutEdge(ii),2};
                [~,edgeNumber] = find(vIDs(1) == IENall);
                
                n = nIntergerPoints(p,0);
                [xi, w] = GaussQuad(n, 1);
                qPoints = BC.getBCQuadrature(xi,w,edgeNumber);
                simplexsf = SimplexShapeFunc(qPoints,p);
                [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();
                
                nssl = length(IENall);
                n = nIntergerPoints(max(pAll),0);
                
                c1Ele = (1/3)*(c1(IENall(1))+c1(IENall(2))+c1(IENall(3)));
                simplex = SimplexStiffMatrix(ele,IENstruct,IENall,n,qPoints,vertexData,ShapeFunc,divShapeFunc,c1Ele);
                
                [elementK,elementF] = simplex.WeakBCOutlettiffMatrix(nout,BCvalEle);
                
                %Assemble the outlet BC into global matrix
                for i = 1:nssl
                    rowTemp = IENall(i);
                    fGlobal(rowTemp) = elementF(i) + fGlobal(rowTemp);
                    for j = 1:nssl
                        columnTemp = IENall(j);
                        indexKglobal = BC.IENindex(rowTemp, columnTemp, row, column);
                        if indexKglobal == 0
                            row = [row, rowTemp];
                            column = [column, columnTemp];
                            value = [value,elementK(i,j)];
                        else
                            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
                        end
                    end
                end
            end
            
            %kGlobal = sparse(row, column, value);
            %K = full(kGlobal)
            %fGlobal
            
        end
        
        
    end
    methods(Static)
        
        
        function [Index] = IENindex(rowTemp, columnTemp, row, column)
            %from a given row column combination, find the correspond index, store the
            %index as another array
            %KIndex(rowTemp,columnTemp) = Index;
            if isempty(row) ==1
                Index = 0;
            else
                indexRow = row == rowTemp;
                indexColumn = column == columnTemp;
                Index = find(indexRow&indexColumn);
                
                %KIndex(rowTemp,columnTemp) = Index;
                
                if isempty(Index)
                    Index = 0;
                else
                    Index = Index;
                end
            end
        end
    end
end
