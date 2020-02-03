classdef SimplexGeometry < handle
    %calculate element level matrix, stiffness matrix and force matrix
    
    properties
        %element ID number
        %EleID;
        %IEN structure including vertexIEN, edgeIEN, faceIEN
        IEN_map;
        %vertex coordinates
        VertexData;
        
    end
    
    methods (Access = public)
        
        function simplex = SimplexGeometry(ele_mesh_ien,vertexCord)
            %contruction
            %simplex.EleID = elementID;
            simplex.IEN_map = ele_mesh_ien;
            simplex.VertexData = vertexCord;
            
        end
        
        function [ JInverse, detJ,gij,x,y,A,hA] = jacobian(tri)
            %calculate jacobian inverse and det J in lambda space, can put
            %this function into private access
            %elementID = tri.EleID;
            vIDs = tri.IEN_map;
            vertexData = tri.VertexData;
            
            %vIDs = IENmesh{elementID,2};
            x1 = vertexData{vIDs(1),2}(1);
            y1 = vertexData{vIDs(1),2}(2);
            x2 = vertexData{vIDs(2),2}(1);
            y2 = vertexData{vIDs(2),2}(2);
            x3 = vertexData{vIDs(3),2}(1);
            y3 = vertexData{vIDs(3),2}(2);
            
            %x1 = 0.2; x2 = 0.4; x3 = 0.2;
            %y1 = 0.6; y2 = 0.6; y3 = 0.8;
            
            x = [x1,x2,x3];
            y = [y1,y2,y3];
            A = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
            
            JInverse = (0.5/A)*[(y3-y1),(x1-x3);(y1-y2),(x2-x1);];
            
            gradx = [(-x1+x2), (-x1+x3);(-y1+y2), (-y1+y3)];
            detJ = 0.5*(gradx(1,1)*gradx(2,2)-gradx(1,2)*gradx(2,1));
            
            dxi1(1) = 0.5/A*(y2-y3); %row
            dxi1(2) = 0.5/A*(x3-x2);
            dxi2(1) = 0.5/A*(y3-y1);
            dxi2(2) = 0.5/A*(x1-x3);
            
            [gij,hA] = tri.metricTensor(dxi1,dxi2,x,y);

        end 
        
    end
    methods(Static)
        function [gij,hA] = metricTensor(dxi1,dxi2,x,y)
            global a;
%             dxidx(1,:) = dxi1;
%             dxidx(2,:) = dxi2;
%             
%             c1=1;%2/sqrt(3);
%             c2=0.5;%1/sqrt(3); % = c1/2
%             
%             K2 = [c1,c2;c2,c1];
%             gij = dxidx'*(K2*dxidx)
%             hA = 1/sqrt((a*(gij*a')));
            
            dNl1 = [-1.0, 1.0, 0.0]; % derivative of N1, N2 and N3 with lambda2
            dNl2 = [-1.0, 0.0, 1.0]; % derivative of N1, N2 and N3 with lambda3
            
            dxdl = zeros(2,2);
            
            for i=1:3
                dxdl(1,1) = dxdl(1,1)+x(i)*dNl1(i); % derivative of x with lambda2
                dxdl(1,2) = dxdl(1,2)+x(i)*dNl2(i); % derivative of x with lambda3
                dxdl(2,1) = dxdl(2,1)+y(i)*dNl1(i); % derivative of y with lambda2
                dxdl(2,2) = dxdl(2,2)+ y(i)*dNl2(i); % derivative of y with lambda3
            end
            
            dldx = inv(dxdl);
            
            c1=1;%2/sqrt(3);
            c2=0.5;%1/sqrt(3); % = c1/2
            
            K2 = [c1,c2;c2,c1];
            gij = dldx'*(K2*dldx);
            hA = 1/sqrt((a*(gij*a')));
            
        end
        
        function [he] = hCalculate(gij,Area)
            global c2;
            global gradMesh;
            global a;
            if (gradMesh == 1)
                if (c2==1)
                    he =2/sqrt(a(1)*gij(1,1)*a(1)+a(1)*gij(1,2)*a(2)+...
                        a(2)*gij(2,1)*a(1)+a(2)*gij(2,2)*a(2));
                elseif (c2==2)
                    he = 2/sqrt(tr(gij)/2);
                end
            else
                he = Area^0.5;
                
            end
        end
            
    
        function [edgeNumber,h] = getBCedge(xCord,yCord,alpha)
            if (abs(xCord(2)-xCord(3))<alpha)
                edgeNumber = 1;
                h = abs(yCord(2)-yCord(3));
            elseif abs(yCord(2)-yCord(3))<alpha
                edgeNumber = 1;
                h = abs(yCord(2)-yCord(3));
            elseif (abs(xCord(1)-xCord(3))<alpha)
                edgeNumber = 2;
                h = abs(yCord(1)-yCord(3));
            elseif abs(yCord(1)-yCord(3))<alpha
                edgeNumber = 2;
                h = abs(xCord(1)-xCord(3));
            elseif (abs(xCord(1)-xCord(2))<alpha )
                edgeNumber = 3;
                h = abs(yCord(1)-yCord(2));
            elseif abs(yCord(1)-yCord(2))<alpha
                edgeNumber = 3;
                h = abs(xCord(1)-xCord(2));
            end
        end
        
        function qPoints = getBCQuadrature(xi,w,edgeNumber)
            if (edgeNumber ==1)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                    qPoints(i,:)=[xi(i),1-xi(i),w(i)];
                end
            elseif (edgeNumber ==2)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                    qPoints(i,:)=[0,1-xi(i),w(i)];
                end
            elseif (edgeNumber ==3)
                lxi = length(xi);
                qPoints = zeros(lxi,3);
                for i = 1:lxi
                    qPoints(i,:)=[1-xi(i),0,w(i)];
                end
            end
        end
    end

       
end