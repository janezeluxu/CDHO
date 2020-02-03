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
        uele;
        utele;
        tauM;
        taunu;
        time;
        
    end
    
    methods (Access = public)
        
        function simplex = quadStiffMatrix(elementID,IEN_mesh,IENall,pAll,nInterger,...
                qPoints,vertexCord,ShapeFunction,divSF,tauEle,variable,variablet,tau1,tau2,timestep)
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
            simplex.uele = variable;
            simplex.utele = variablet;
            simplex.tauM = tau1;
            simplex.taunu = tau2;
            simplex.time = timestep;
            %simplex.SF1 = tauShapeFunc;
            
        end
        
        function [elementK,elementF,tauuse] = eleStiffMatrix(simplex,JInv, detJall,gijall,itau,cc,quadraturePoints,DC)
            %element stiffness matrix and force
            
            %global stabflag; 
            global kappa;
            global a; 
            global casenumber; 
            %global IntByPart;
            
            global dt;
            %global gamma;
            %global af;
            global am;
            %cc = (af*gamma*dt);
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
            u_ele = simplex.uele;
            ut_ele = simplex.utele;
            %taumquad = simplex.tauM;
            %taunuquad = simplex.taunu ;
            t = simplex.time;  
            %quadraturePoints = QuadGaussPoints(n);
            nP = size(quadraturePoints,1);
            sizeN = length(IENall);
            elementK = zeros(sizeN, sizeN);
            elementF = zeros(sizeN,1);
            
            vIDs = ien_mesh;
            x1 = vertexData{vIDs(1),2}(1);
            y1 = vertexData{vIDs(1),2}(2);
            x2 = vertexData{vIDs(2),2}(1);
            y2 = vertexData{vIDs(2),2}(2);
            x3 = vertexData{vIDs(3),2}(1);
            y3 = vertexData{vIDs(3),2}(2);
            x4 = vertexData{vIDs(4),2}(1);
            y4 = vertexData{vIDs(4),2}(2);
            
            x = [x1,x2,x3,x4];
            y = [y1,y2,y3,y4];
            
            hA = sqrt((x1-x2)^2+(y1-y2)^2);
            A = hA^2;

            count = 1;
            tauquad = zeros(nP,1);
            for k = 1:nP 
                JInverse = reshape(JInv(count:count+3),2,2);
                detJ = detJall(k);
                gij = reshape(gijall(count:count+3),2,2);
                count = count+4;
                
                xi = quadraturePoints(k,1);
                eta = quadraturePoints(k,2);
                
                %[JInverse, detJ,gij,x,y,A,hA] =  ...
                %    getjacobian(ien_mesh,vertexData,xi,eta);
                Jw = detJ*quadraturePoints(k,3);
                
                N1 = (1-xi)*(1-eta);
                N2 = xi*(1-eta);
                N4 = xi*eta;
                N3 = (1-xi)*eta;
                
                xc = N1*x(1)+N2*x(2)+N3*x(4)+N4*x(3);
                yc = N1*y(1)+N2*y(2)+N3*y(4)+N4*y(3);
                
                sourceTerm = MMS(t,xc,yc,casenumber);
                
                if isa(a{1}, 'function_handle')==1
                    velocity1 = a{1}(xc,yc);
                else
                    velocity1 = a{1};
                end
                if isa(a{2}, 'function_handle')==1
                    velocity2 = a{2}(xc,yc);
                else
                    velocity2 = a{2};
                end
                
                velocity = [velocity1,velocity2];
                
                gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
                NbGlobal = ShapeFunc(:,k);
                SF = ShapeFunc(:,k);
                SFT = SF';
                
                %uquad = SFT*u_ele;
                %u_ele
                diffu = gradNaGlobal'*u_ele;
                
                ut = SFT*ut_ele;
                ru = ut+velocity*diffu-sourceTerm;
                
                %tauMquad = tauMele;
                %tau = simplex.tauFunc(gij,hA,Plist,taumquad,velocity);
                %tau = taumquad(k);
                %tau2 = simplex.tau2Func(tau,ru,gij,SF,diffu,Plist);
                if(itau == 2)
                    p = min(Plist);
                    gij_u = 4*gij*p^2;
                    tau1sqinv = velocity*gij_u*velocity';
                    static_factor = (3)^2;
                    tau2sqinv = static_factor*sum(dot(kappa*gij_u,kappa*gij_u));
                    C1 = 1;
                    taut = (2*C1/dt)^2;
                    tau = 1/sqrt((taut+tau1sqinv+tau2sqinv));
                elseif(itau == 3)
                    tau = tauMele;
                    
                elseif(itau == 0)
                    tau = 0;
                end
                
                tauquad(k) = tau;
                if (DC == true)
                    gij = 4*gij;
                    Invgij = inv(gij);
                    denominator = (diffu'*Invgij*diffu);                    
                    if denominator>0
                        miu1 = sqrt(ru^2/denominator);
                        miu2 = tau*ru^2/denominator;
                        miuh = max(0,miu1-miu2);
                        miu = miuh;
                    else
                        miu = 0;
                    end
                    tau2 = miu*(Invgij);
                else
                    tau2 = 0;
                end
                
                elementF = elementF - ShapeFunc(:,k)*Jw*sourceTerm;
                elementF = elementF + ShapeFunc(:,k)*ut*Jw;
                elementF = elementF + ShapeFunc(:,k)*velocity*diffu*Jw;
                elementF = elementF + gradNaGlobal*kappa*diffu*Jw;
                elementF = elementF + gradNaGlobal*velocity'*tau*ru*Jw;
                elementF = elementF + gradNaGlobal*tau2*diffu*Jw;
                
                elementK = elementK + SF*SFT*am*Jw+gradNaGlobal*kappa*gradNaGlobal'*Jw*cc;
                % Advection
                %if (strcmp(IntByPart,'Yes')==1)
                %    elementK = elementK-gradNaGlobal*velocity'*NbGlobal'*Jw*cc;
                %else
                elementK = elementK+NbGlobal*(velocity*gradNaGlobal')*Jw*cc;
                %end
                % Stabilizer
                %if(stabflag ==1)
                elementK = elementK + (gradNaGlobal*velocity')*tau*...
                    (NbGlobal'*am+velocity*gradNaGlobal'*cc)*Jw;
                elementK = elementK + gradNaGlobal*tau2*gradNaGlobal'*Jw*cc;
                %end
            end
            elementF = -elementF;
            tauuse = max(tauquad);
        end

        function [ tau,taut,tau1sqinv,tau2sqinv] = tauFunc(simplex,gij,h,Plist,tauMquad,a)
            global itau;
            global nsd; 
            global kappa;
            global dt; 
            if(itau == 2)
                p = min(Plist);
                gij = 4*gij*p^2;
                tau1sqinv = a*gij*a';
                static_factor = (3)^2;
                tau2sqinv = static_factor*sum(dot(kappa*gij,kappa*gij));
                C1 = 1;
                taut = (2*C1/dt)^2;
                
                tau = 1/sqrt((taut+tau1sqinv+tau2sqinv));
                
            elseif(itau == 3)
                tau = tauMquad;
                
            elseif(itau == 0)
                tau = 0;
            end
            
        end
        
        function [ tau] = tau2Func(simplex,tau,ru,gij,SF,diffu,Plist)  
            Invgij = inv(gij);
            denominator = (diffu'*Invgij*diffu);          
            miu1 = sqrt(ru^2/denominator);
            miu2 = tau*ru^2/denominator;
            
            miu3 = 2*tau*ru^2/denominator;
            miuh = max(0,miu1-miu2);
            miu = min(miuh,miu3);
            miu = 0;
            tau = miu*(Invgij);
        end
        
    end
    
    methods(Static)
        
        function [Source] = MMS(t,x,y,casenumber)
            if casenumber == 300
                x = 10*x-1;
                y = 10*y-1;
                if (sqrt(x^2+y^2)>1)
                    Source = 0;
                else
                    Source = exp(-t^1)*cos((pi/2)*sqrt(x^2+y^2));
                end
            elseif casenumber == 400
                Source = 0;
            end
            
        end
    end
end