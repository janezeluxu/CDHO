classdef geo_Mesh < handle
    
    properties (Access = private)
        %number of Points in x direction
        meshien;
        nodecoord;
        BCD;
        
        % Order p list
        Order;
        PV;
        PH;
        PT

    end
    
    methods(Access = public)
        
        function geoMesh = geo_Mesh(meshIEN,...
                vertexCood,BC_array,varargin)  
        
        geoMesh.meshien = meshIEN;
        geoMesh.nodecoord = vertexCood;
        geoMesh.BCD = BC_array;
        
        if length(varargin) ==1
            geoMesh.Order = varargin{1}; 
            geoMesh.PV = nan;
            geoMesh.PH = nan;
            geoMesh.PT = nan;
        elseif length(varargin) ==3
            geoMesh.PV = varargin{1};
            geoMesh.PH = varargin{2};
            geoMesh.PT = varargin{3};
            geoMesh.Order = nan;
        else
          disp('number of input is wrong!')
        end
        end
        
        
        function [meshData] = MeshData(geoMesh)
            %[faceID, vertexID]
            global Grid_size;
            meshIENarray = geoMesh.meshien;
            meshData = cell(Grid_size,2);
            %element number and vertex
            for i  = 1:Grid_size
                meshData{i,1} = i; 
                meshData{i,2} = meshIENarray(i,:);
                %meshData{i,2}(1) = meshIENarray(i,1);
                %meshData{i,2}(2) = meshIENarray(i,3);
                %meshData{i,2}(3) = meshIENarray(i,2);
            end
        end
        
        function [vertexData] = vertexMesh(geoMesh,meshData)
            %[vID, vCoordinate, element around vertex, vertex around vertex]
            global node_Size;
            vertexCood = geoMesh.nodecoord;
            vertexData = cell(node_Size,2);
            for vID  = 1:node_Size
                vertexData{vID,1} = vID;
                vertexData{vID,2} = vertexCood(vID,:);
            end
            
            vertexFace = cell2mat(meshData(:,2));
            for i = 1:length(vertexData)
                [faceID] = geoMesh.FaceVertex(i,vertexFace);
                vertexData{i,3} = unique(faceID'); %vertex face relation
                vvID = [];
                for j = 1:length(faceID)
                    vvID = [vvID,meshData{faceID(j),2}];
                end
                vertexData{i,4} = unique(vvID);
            end            
        end
        
        function [Area] = elementArea(geoMesh,meshData,vertexData)
            global Grid_size;
            Area = zeros(Grid_size,1);
            for i = 1:Grid_size
                [vID] = meshData{i,2};
                vCood1 = vertexData{vID(1),2};
                x1 = vCood1(1);
                y1 = vCood1(2);
                vCood2 = vertexData{vID(2),2};
                x2 = vCood2(1);
                y2 = vCood2(2);
                vCood3 = vertexData{vID(3),2};
                x3 = vCood3(1);
                y3 = vCood3(2);
                vCood4 = vertexData{vID(4),2};
                x4 = vCood4(1);
                y4 = vCood4(2);
                hh = (x2-x1)^2+(y2-y1)^2;
                %Area(i) = (x2-x1)*(y3-y1);
                Area(i) =hh;
            end
        end
            
        function [IBC,BCval] = BoundaryCondition(geoMesh)
            %% All boundary DOFs list assuming all BCvalue are constant not
            %%functions
            global Dirichlet_BC;
            
            IBC = [];
            BCval = [];
            IDall = geoMesh.BCD;
            for i = 1:size(IDall,1)
                vID = i;
                v_tag = IDall(i);
                for j = 1:size(Dirichlet_BC{1},1)
                    if (v_tag == Dirichlet_BC{1}(j))
                        IBC = [IBC,vID];
                        BCval = [BCval,Dirichlet_BC{2}{j}];
                    end
                end
            end
        end
        
        
        
    end

    methods(Static)
        
        function [edgeVertexList] = edgeVertex(edgeID, meshData,meshDataEdge)
            %from the edgeID, find the two vertexID
            [rowID,columnID] = find((edgeID==meshDataEdge),1);
            edgeVertexList = [meshData{rowID,2}(columnID),meshData{rowID,2}(mod(columnID,3)+1)];
            
        end
        
        function [faceID] = FaceEdge(edgeID,edgeFace)
            %from edgeID, find the two/one face 
            [faceID,~] = find(edgeID == edgeFace);
            %faceID = row;
        end
        
        function [faceID] = FaceVertex(vID,vertexFace)
            %from vID, find the face 
            [faceID,~] = find(vID == vertexFace);
            %faceID = row;
        end
        
        function [BCval] = AssignBCval(BCI,vertexID,edgeDOF,BCval,vertexData)
            global BC;
            global StepFunc;

            if (strcmp(BC(BCI),'StepFunc')==1)
                disp ('StepFunc!!!');
                count = 1;
                for i = 1:(length(vertexID))
                    vcord = vertexData{vertexID(i),2};
                    
                    if (BCI == 3 || BCI == 4)
                        
                        if vcord(2)<StepFunc(3) %stepfunction in y direction
                            BCval(count) = StepFunc(1);
                        else
                            BCval(count) = StepFunc(2);
                        end
                    else
                        if vcord(1)<StepFunc(3)%stepfunction in x direction
                            BCval(count) = StepFunc(1);
                        else
                            BCval(count) = StepFunc(2);
                        end
                    end
                    count = count+1;
                end
            else            
            count = 1;
            for i = 1:(length(vertexID))
                BCval(count) = BC{BCI};
                count = count+1;
            end
            for i = 1:length(edgeDOF)
                BCval(count) = BC{BCI};
                count = count+1;
            end
            end
        end
        
        
    end
end