classdef IEN < handle
    
    properties (Access = private)
        %input from geo_Mesh class
        MeshData; 
        solutionData; 
        orderData;
        usage;
        
    end
    
    methods(Access = public)
        function ien = IEN(mesh_ien,solution_ien,order_list,face_edge)
            ien.MeshData = mesh_ien;
            ien.solutionData = solution_ien;
            ien.orderData = order_list;
            ien.usage = face_edge;
            
        end
        function [solutionIEN,pOrder] = Construct_IEN(ien)
            solutionIEN = ien.solutionData;
            pOrder = ien.orderData;
            face_usage = ien.usage;
            %vertexIEN[face, vertexDOF]
            %edgeIEN [edge1Order, edge1DOF],[edge2Order,
            %edge2DOF],[edge3Order, edge3DOF]
            %faceIEN [faceID] [faceOrder, faceDOF]
%             global Grid_size
%             Grid_size
%             meshData = ien.MeshData;
%             solution = ien.solutionData;
%             order = ien.orderData;
%             
%             vertexIEN = cell(Grid_size,3);
%             edgeIEN = cell(Grid_size,3);
%             faceIEN = cell(Grid_size,3);
%             
%             %loop over all elements, find the correspond vID, edgeID, faceID
%             for i  = 1:Grid_size
%                 vertexIEN{i,1} = [meshData{i,1}];
%                 vertexIEN{i,2} = [meshData{i,2}];
%                 
%                 %% !!!face uses edge in reverse order, will add usage information later
%                  edgeIEN{i,1} = [order(i,1), solution(i,4:order(i,1)+3)];
%                  edgeIEN{i,2} = [order(i,2),solution(i,:)];
%                  edgeIEN{i,3} = [order(i,3), solution(i,:)];
%                               
%                 fID = i;
%                 faceIEN{i,1} = [meshData{i,1}];
%                 faceIEN{i,2} = [order(i,4),solution(i,:)];
%                 
%             end
%             %IENstruct = struct('vertexIEN',vertexIEN,'edgeIEN',edgeIEN,'faceIEN',faceIEN)
%             solution(1,4:order(1,1)+2)
        end

    end
    
    
end