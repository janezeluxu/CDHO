classdef quadStiffMatrix < handle
    %calculate element level matrix, stiffness matrix and force matrix
    
    properties
        %element ID number
        EleID;
        %element level IEN arra
        IENs;
        orderList;
        mesh_IEN;
        %number of intergral points;
        nInt;
        %vertex coordinates
        VertexData;
        %element shape functions
        EleShapeFunc;
        %element shape function derivitives
        EledivShapeFunc;
        tauElement;
        Quadrature;
        
    end
    
    methods (Access = public)
        
        function simplex = quadStiffMatrix(elementID,IEN_mesh,IENall,pAll,nInterger,...
                qPoints,vertexCord,ShapeFunction,divSF,tauEle)
            %contruction
            simplex.EleID = elementID;
            simplex.mesh_IEN = IEN_mesh;
            simplex.IENs = IENall;
            simplex.orderList = pAll;
            simplex.nInt = nInterger;
            simplex.VertexData = vertexCord;
            simplex.EleShapeFunc = ShapeFunction;
            simplex.EledivShapeFunc = divSF;
            simplex.tauElement = tauEle;
            simplex.Quadrature = qPoints;
            %simplex.SF1 = tauShapeFunc;
            
        end
        
        function [elementK,elementF] = eleStiffMatrix(simplex)
            %element stiffness matrix and force
            
            global stabflag; 
            global kappa;
            global a; 
            global force; 
            global IntByPart;
            
            IENall = simplex.IENs;
            Plist = simplex.orderList;
            ShapeFunc = simplex.EleShapeFunc;
            divShapeFunc = simplex.EledivShapeFunc;            
            n = simplex.nInt;
            eleID = simplex.EleID;
            ien_mesh = simplex.mesh_IEN;
            vertexData = simplex.VertexData;
            tauMele = simplex.tauElement;
            %ShapeFunctau = simplex.SF1;
            
            %quadraturePoints = TriGaussPoints(n);   
            quadraturePoints = QuadGaussPoints(n);
            nP = size(quadraturePoints,1);
            sizeN = length(IENall);
            elementK = zeros(sizeN, sizeN);
            elementF = zeros(sizeN,1);
            
            sourceTerm = zeros(nP,1);
            %sourceTerm = kappa*simplex.MMS_Diffusion(x,y,force)+a*simplex.MMS_Adv(x,y,force);
            
            
            for k = 1:nP  
                xi = quadraturePoints(k,1);
                eta = quadraturePoints(k,2);
                
                [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                    getjacobian(ien_mesh,vertexData,xi,eta);
                Jw = detJ*quadraturePoints(k,3);
                
                elementF = elementF+ShapeFunc(:,k)*Jw*sourceTerm(k);
                gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
                
                % Stability part of force vector
                if(stabflag==1)
                    % tau calculation
                    %fileID = fopen('./testvtk/taustatic.txt','a+');
                    %SFtau = ShapeFunctau(:,k);
                    %tauMquad = SFtau'*tauMele;
                    tauMquad = tauMele;
                    tau = simplex.tauFunc(gij,hA,Plist,tauMquad);
                    %fprintf(fileID,'%2.16f \n',tau);
                    %fclose(fileID);
                    elementF = elementF + Jw*gradNaGlobal*a'*tau*sourceTerm(k);
                end
                
                NbGlobal = ShapeFunc(:,k);
                
                %diffusion
                elementK = elementK + gradNaGlobal*kappa*gradNaGlobal'*Jw;
                % Advection
                if (strcmp(IntByPart,'Yes')==1)
                    elementK = elementK-gradNaGlobal*a'*NbGlobal'*Jw;
                else 
                    elementK = elementK+NbGlobal*(a*gradNaGlobal')*Jw;
                end
                % Stabilizer                
                if(stabflag ==1)
                    elementK = elementK + (gradNaGlobal*a')*tau*(gradNaGlobal*a')'*Jw;
                end
            end
            
        end

        function [ tau] = tauFunc(simplex,gij,h,Plist,tauMquad)
            global itau;
            global nsd; 
            global kappa;
            global a; 
            if(itau == 2)
                p = min(Plist);
                gij = 4*gij*p^2;
                tau1sqinv = a*gij*a';
                static_factor = (3)^2;
                tau2sqinv = static_factor*sum(dot(kappa*gij,kappa*gij));
                tau = 1/sqrt((tau1sqinv+tau2sqinv));
                
            elseif(itau == 3)
                tau = tauMquad;
                
            elseif(itau == 0)
                tau = 0;
            end
            
        end
        
    end
end