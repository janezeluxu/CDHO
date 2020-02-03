classdef SimplexStiffMatrix < handle
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
        
        function simplex = SimplexStiffMatrix(elementID,IEN_mesh,IENall,pAll,nInterger,...
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
            
            quadraturePoints = TriGaussPoints(n);            
            nP = size(quadraturePoints,1);
            sizeN = length(IENall);
            elementK = zeros(sizeN, sizeN);
            elementF = zeros(sizeN,1);
            
            %divShapeFunc;
            tri = SimplexGeometry(ien_mesh,vertexData);
            [JInverse, detJ,gij,xCord,yCord,A,hA] =  tri.jacobian();
            %if(eleID==0)
            %    JInverse, detJ,gij,xCord,yCord,A,hA
            %end
            %h = tri.hCalculate(gij,A);
            
            Jw = detJ*quadraturePoints(:,3);
            
            lam = [1-quadraturePoints(:,1)-quadraturePoints(:,2),quadraturePoints(:,1:2)];
            x = xCord*lam';
            y = yCord*lam';
            
            sourceTerm = simplex.MMS(x,y,force);
            %sourceTerm = kappa*simplex.MMS_Diffusion(x,y,force)+a*simplex.MMS_Adv(x,y,force);
            
            for k = 1:nP  
                elementF = elementF+ShapeFunc(:,k)*Jw(k)*sourceTerm(k);
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
                    elementF = elementF + Jw(k)*gradNaGlobal*a'*tau*sourceTerm(k);
                end
                
                NbGlobal = ShapeFunc(:,k);
                
                %diffusion
                elementK = elementK + gradNaGlobal*kappa*gradNaGlobal'*Jw(k);
                % Advection
                if (strcmp(IntByPart,'Yes')==1)
                    elementK = elementK-gradNaGlobal*a'*NbGlobal'*Jw(k);
                else 
                    elementK = elementK+NbGlobal*(a*gradNaGlobal')*Jw(k);
                end
                % Stabilizer                
                if(stabflag ==1)
                    elementK = elementK + (gradNaGlobal*a')*tau*(gradNaGlobal*a')'*Jw(k);
                end
            end
            
        end

        function [ tau] = tauFunc(simplex,gij,h,Plist,tauMquad)
            global itau;
            global nsd; 
            global kappa;
            global a; 
            %global static_factor;
            %h = Area^0.5;
            if(itau == 1)
                alpha = norm(a)*h/(2*norm(kappa));
                zi = coth(alpha)-1/alpha;
                if(nsd==1)
                    tau = h*zi/(2*norm(a));
                else
                    %this does not work
                    tau = h*zi*norm(a)/norm(kappa);
                end
            elseif(itau == 2)
                if(nsd==1)
                        tau1 = h/(2*a(1));
                        tau2 = h^2/(4*kappa(1,1));
                        tau = 1/sqrt(tau1^-2 + 9*tau2^-2);
                else
                    gij = 4*gij;
                    tau1sqinv = a*gij*a';
                    p = min(Plist);
                    static_factor = (3*min(Plist)^2)^2;
                    tau2sqinv = static_factor*sum(dot(kappa*gij,kappa*gij));
                    tau = 1/sqrt((tau1sqinv+tau2sqinv));
                    
                end
                
            elseif(itau == 3)
                tau = tauMquad;
                
            elseif(itau == 0)
                tau = 0;
            end
            
        end
        
    end
    methods(Static)
        
        function [Source] = MMS(x,y,q)
            global kappa;
            global a;
            global direction;
            
            if (strcmp(direction,'Y')==1)
                x = y;
            end              
            if q ==0
                Source = zeros(length(x));
            elseif q == 1
                Source = 1*ones(length(x));                
            elseif q == 2
                Source_diff = kappa(1,1)*(-2*(x.^2-1))+kappa(2,2)*(-2*(y.^2-1));
                Source_adv = a(1)*2*x.*(y.^2-1)+a(2)*2*y.*(x.^2-1);
                Source = Source_diff+Source_adv;
            elseif q ==3
                Source = 12*x.^2;
            elseif q ==4
                Source = 2*sin(x).*sin(y);
                
            elseif q == 5
                Source = zeros(1,length(x));
                for i = 1:length(x)
                        if x(i)<= -0.6                           
                            Source = 0*ones(length(x),1);
                        elseif x(i) >0.6
                            Source = 0*ones(length(x),1);
                        else
                            Source = (5/3)*ones(length(x),1);
                        end
                end
                
                elseif q == 6
                Source = zeros(1,length(x));
                for i = 1:length(x)
                        if x(i)<= -0.6
                            Source = 0*ones(length(x),1);
                        elseif x(i) >0.6
                            Source = 0*ones(length(x),1);
                        else
                            ax = 1/(4*0.6^3);
                            Source(i) = 12*ax*x(i)^2;
                        end
                end
            elseif q == 7
                Source = zeros(1,length(x));
                for i = 1:length(x)
                    if x(i)<= -0.6
                        Source(i) = 2*(x(i)+1);
                    elseif x(i) >0.6
                        Source(i) = -2*(x(i)-1);
                    else
                        Source(i) = 2*((-5/6)*x(i)^2+0.7)+(5/3)*(1-y(i)^2);
                    end
                end
            
            elseif q == 8
                Source = zeros(1,length(x));
                for i = 1:length(x)
                    if (x(i)> -0.6 && x(i)< 0.6 && y(i)> -0.6 && y(i)< 0.6)
                        Source(i) = -2*(x(i)^2+y(i)^2-0.72);
                    else
                        Source(i) = 0;
                    end
                end
            end
        end

        
    end
       
end

