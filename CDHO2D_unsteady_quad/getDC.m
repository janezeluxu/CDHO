function [taus,tau_DC,taut_percent,tauadv_percent,taudiff_percent] = ...
    getDC(t,variable,variable_t,vertexData,JInv,detJall,gijall,meshData,solution_ien,p_ien,tauM,p)
global Grid_size;
global a;
global casenumber;
global itau;
global kappa;
global dt;
global DC;
%n = nIntergerPoints(p,0);
%qPoints = TriGaussPoints(n);
n = nIntergerPoints(p,0);
qPoints = QuadGaussPoints(n);

nP = size(qPoints,1);
tau_DC = zeros(Grid_size,1);
taus = zeros(Grid_size,1);
taut_percent = zeros(Grid_size,1);
tauadv_percent = zeros(Grid_size,1);
taudiff_percent = zeros(Grid_size,1);
%[ShapeFuncTable,divSFtable] = ShapeTable(0,p);
[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);
for ele = 1:Grid_size
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    nssl = length(IENall);
    tauele = zeros(4,1);
    for i = 1:4
       tauele(i) = tauM(IEN_mesh(i));
    end
    tauMele = max(tauele);
    
    u_ele = zeros(nssl,1);
    ut_ele = zeros(nssl,1);
    for i = 1:nssl
        u_ele(i) = variable(IENall(i));
        ut_ele(i) = variable_t(IENall(i));
    end
    taut_p = zeros(nP,1);
    tauadv_p = zeros(nP,1);
    taudiff_p = zeros(nP,1);
    tauDC = zeros(nP,1);
    tau1 = zeros(nP,1);
    
    vIDs = IEN_mesh;
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
       
    JInverseele = JInv(ele,:);
    detJele = detJall(ele,:);
    gijele = gijall(ele,:);
    
    count = 1;
    for k = 1:nP
        
        JInverse = reshape(JInverseele(count:count+3),2,2);
        detJ = detJele(k);
        gij = reshape(gijele(count:count+3),2,2);
        count = count+4;
        
        xi = qPoints(k,1);
        eta = qPoints(k,2);
        
        %Jw = detJ*quadraturePoints(k,3);
        
        N1 = (1-xi)*(1-eta);
        N2 = xi*(1-eta);
        N4 = xi*eta;
        N3 = (1-xi)*eta;
        
        xc = N1*x(1)+N2*x(2)+N3*x(4)+N4*x(3);
        yc = N1*y(1)+N2*y(2)+N3*y(4)+N4*y(3);
        
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
        
        sourceTerm = MMS(t,xc,yc,casenumber);
        
        SF = ShapeFunc(:,k);
        SFT = SF';
        gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
        diffu = gradNaGlobal'*u_ele;
        ut = SFT*ut_ele;
        ru = ut+velocity*diffu-sourceTerm;
        
        tauMquad = tauMele;
        if(itau == 2)
            order = min(pAll);
            gij_u = 4*gij*order^2;
            tau1sqinv = velocity*gij_u*velocity';
            static_factor = (3)^2;
            tau2sqinv = static_factor*sum(dot(kappa*gij_u,kappa*gij_u));
            C1 = 1;
            taut = (2*C1/dt)^2;
            tau = 1/sqrt((taut+tau1sqinv+tau2sqinv));
            
            taut_p(k) = taut/(taut+tau1sqinv+tau2sqinv);
            tauadv_p(k) = tau1sqinv/(taut+tau1sqinv+tau2sqinv);
            taudiff_p(k) = tau2sqinv/(taut+tau1sqinv+tau2sqinv);
        elseif(itau == 3)
            tau = tauMquad;
            
        elseif(itau == 0)
            tau = 0;
        end
        
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
            D = eig(tau2,'matrix');
            %miu_12 = D(1,2);
            miu_22 = D(2,2);
            
            tauDC(k) = miu_22;
        
        else
            tauDC(k) = 0;
        end
        tau1(k) = tau;
    end
    taut_percent(ele) = max(taut_p);
    tauadv_percent(ele) = max(tauadv_p);
    taudiff_percent(ele) = max(taudiff_p);
    tau_DC(ele) = max(tauDC);
    taus(ele) = max(tau1);
end
end

function [ tau] = tauFunc(gij,Plist,tauMquad)
global itau;
global kappa;
global a;
%global static_factor;
%h = Area^0.5;
if(itau == 2)
%     gij = 4*gij;
%     tau1sqinv = a*gij*a';
%     p = min(Plist);
%     static_factor = (3*min(Plist)^2)^2;
%     tau2sqinv = static_factor*sum(dot(kappa*gij,kappa*gij));
%     tau = 1/sqrt((tau1sqinv+tau2sqinv));

    gij = 4*gij;
    p = min(Plist);
    tau1sqinv = a*gij*a'*p^2;
    
    static_factor = (3)^2;
    tau2sqinv = static_factor*sum(dot(kappa*gij*p^2,kappa*gij*p^2));
    tau = 1/sqrt((tau1sqinv+tau2sqinv));
    
    
elseif(itau == 3)
    tau = tauMquad;
end


end
        
function [ miu] = tau2Func(tau,ru,gij,diffu,pAll)
Invgij = inv(gij);
denominator = (diffu'*Invgij*diffu);
p = min(pAll);
if denominator>0
    miu1 = (1/p)*sqrt(ru^2/denominator);
    miu2 = tau*ru^2/denominator;
    miuh = max(0,miu1-miu2);
    
    %miu3 = 2*miu2;
    %miu = min(miuh,miu3);
    miu = miuh;
else
    miu = 0;
end
%miu = 0;
%tauquad = miu*(Invgij);
end