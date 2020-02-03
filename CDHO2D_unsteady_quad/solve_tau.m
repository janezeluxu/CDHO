function [tau] = solve_tau(t,u_n,ut_n,tau,vertexData,...
                meshData,solution_ien,order_ien, ...
                JInverse, detJ,gijG,order,IBC,Area,...
                patchArea,patchOrder,uHBCnode,istp,aref)   
global node_Size;
global Grid_size;
global itau;
global casenumber;
%% solve tau_m
[tau_g,c1node,c2node,L,M] = solvetauHOmap(t,order,IBC,vertexData,meshData,...
    solution_ien,order_ien,JInverse, detJ,...
    u_n,ut_n,patchArea,patchOrder,uHBCnode,Area,aref);
%tau_g
%% update taum element based
Tavg =  getTave_ele(vertexData,meshData,solution_ien,order_ien,gijG);
%Tavg = ones(Grid_size,1)*1e-5;
tau = relxation_ele(tau,tau_g,Tavg);

% 
% % %% update taum solving scalar relaxation problem
% % IBCI = IBC_tau{1};
% % updateTm = updateTau(t,Hhfactor,L,M,Tavg,variable,vertexData,IBCI,...
% %     order,Area,rowI,columnI,KindexI,rowIndexI,columnIndexI,...
% %     IENmesh,IENall,pAll,SF, gradNG,JwG,gijG,sizea,sizeb,sizec,sized,sizee,sizef);
% % 
% % [Ilm_np1,Itlm_np1, Imm_np1,Itmm_np1] = getI(updateTm,...
% %     Ilm_n,Itlm_n,Ilm_np1,Itlm_np1,Imm_n,Itmm_n,Imm_np1,Itmm_np1);
% % 
% % c1 = updateTm.getTm(Ilm_np1,Imm_np1,IBCI);
% 
% % %% update c1 sudo time step problem
% % delt = 3/7;
% % c1 = (1/(1+delt))*c1+(delt/(1+delt))*c1node_g;
% 
% %% update c1 using relaxation on c1
% Tavg =  getTave(vertexData,gijG);
% %Tavg = ones(node_Size,1);
% tau = relxation(tau,tau_g,Tavg);
% 
% %% update c1 using relaxation on L and M
% [LM,MM] = getLM(L,M);
% Tavg =  getTave(vertexData,IENmesh,IENall,pAll,gijG,variable);
% Ilm_np1 = relxation(Ilm_n,LM,Tavg);
% Imm_np1 = relxation(Imm_n,MM,Tavg);
% IBCI = IBC_tau{1};
% [c1] = updatec1(Ilm_np1,Imm_np1,IBCI);
% 
%% update taum element based
% Tavg =  getTave_ele(vertexData,meshData,solution_ien,order_ien,gijG);
% %Tavg = ones(Grid_size,1)*1e-5;
% tau = relxation_ele(tau,tau_g,Tavg)
% % [maxtausolve_tau,maxIndex] = max(taum);
% % maxtausolve_tau
% % [mintau_solvetau,minIndex] = min(taum);
% % mintau_solvetau
% %% output taum
% output_filename = strcat('./testvtk/cases/cylindertest/',itaucase,'p1/',casenum,string(istp),'tau.vtk')
% writevtkfile(IENmesh,vertexData,[c1node_g(1:node_Size,1),...
%     c2node(1:node_Size)],c1(1:node_Size),zeros(Grid_size,1),zeros(Grid_size,1));
% twod_to_vtk_test (output_filename);
% 

tauELE = zeros(Grid_size,1);
for ele = 1:Grid_size 
    tauELE(ele) = max(tau(ele,:));
end

if mod(istp,10) == 0
    output_filename = strcat('./testvtk/cases/pulse/dynamicp7/dynamic',string(istp),'tau.vtk')
    writevtkfile(meshData,vertexData,L(1:node_Size,1),M(1:node_Size,1),tauELE);
    twod_to_vtk_test (output_filename);
end
end

function [c1_np1] = relxation(c1node_g,c1_n,Tavg)
global dt;
global node_Size;
%Tavg = dt*ones(node_Size);
c1_np1 = zeros(node_Size,1);
for node = 1:node_Size 
    c1_np1(node)= (c1_n(node)-c1node_g(node))*exp(-dt/Tavg(node))+c1node_g(node);
end
end

function [taum_np1] = relxation_ele(taum,taum_g,Tavg)
global dt;
global Grid_size;
%Tavg = dt*ones(node_Size);
taum_np1 = zeros(Grid_size,4);
for ele = 1:Grid_size 
    for i = 1:4
       taum_np1(ele,i) = (taum(ele,i)-taum_g(ele,i))*exp(-dt/Tavg(ele))+taum_g(ele,i);
    end
end
end

function [taumnode] = updatec1(Ilm,Imm,IBCtau)
global node_Size
taumnode = zeros(node_Size,1);

for node = 1:node_Size
    if norm(Imm(node))<1e-10
        taumnode(node) = 0;
    else
        taumnode(node) = Imm(node)\Ilm(node);
    end
end

for node = 1:node_Size
    if taumnode(node) <= 0
        %taumnode(node) = 0;
        taumnode(node)=abs(taumnode(node));
    end
end

for i = 1:length(IBCtau)
    taumnode(IBCtau(i)) = 0;
end

end

function [LM,MM] = getLM(L,M)
global node_Size;

LM = zeros(node_Size,1);
MM = zeros(node_Size,1);
%testLM= zeros(TotalNode,1);
for i = 1:node_Size
    LM(i) = (L(i,1)*M(i,1)+L(i,2)*M(i,2)+L(i,3)*M(i,3)+L(i,4)*M(i,4));
    MM(i) = (M(i,1)*M(i,1)+M(i,2)*M(i,2)+M(i,3)*M(i,3)+M(i,4)*M(i,4));
    %testLM(i) = MM(i)\LM(i);
    if (LM(i)*MM(i)<0)
        LM(i) = 0;
        %LM(i) = abs(LM(i));
        %MM(i) = abs(MM(i));
    end
end
end

function tauM = gettaum(Area,IENmesh,c1,c2)
global Grid_size
tauM = zeros(Grid_size,1);
for ele = 1:Grid_size 
    %% use maximum nodal taum as element taum
    IEN_mesh = IENmesh(ele,:);
    nsf = length(IEN_mesh);
    hh = sqrt(Area(ele));
    taum_ele = zeros(nsf,1);
    for i = 1:nsf
        c2_ele = c2(IEN_mesh(i));
        c1_ele = c1(IEN_mesh(i));
        taum_ele(i) = c1_ele*hh^c2_ele;
    end
    tauM(ele) = max(taum_ele);
end
end

function TT =  getTave(vertexData,gijG)
global node_Size
global dt;
global a;
TT = zeros(node_Size,1);
for node = 1:node_Size  
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);  
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        gijele = gijG(ele,:);
        gij = reshape(gijele(1:4),2,2);
    end
    xc = vertexData{node,2}(1);
    yc = vertexData{node,2}(2);           
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
    
    aa = [velocity1,velocity2];
    
    aT = aa';
    tau1sqinv = aa*gij*aT;
    C1 = 1;
    taut = (2*C1/dt)^2;
    taus = 1/sqrt((taut+tau1sqinv));
    TT(node) = 10*taus;
end
end

function TT_ele =  getTave_ele(vertexData,meshData,solution_ien,p_ien,gijG)
global Grid_size
global a;
global dt;
TT_ele = zeros(Grid_size,1);
for ele = 1:Grid_size
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    
    nssl = length(IENall);
    n = nIntergerPoints(max(pAll),0);
    qPoints = QuadGaussPoints(n);
    
    gij_ele = gijG(ele,:);
    nP = size(qPoints,1);
    TT = zeros(nP,1);
    
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
    count = 1;
    
    for k = 1:nP
        xi = qPoints(k,1);
        eta = qPoints(k,2);
        %% scalar
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
        
        aa = [velocity1,velocity2];
        
        aT = aa';
        gij = reshape(gij_ele(count:count+3),2,2);
        count = count+4;
        tau1sqinv = aa*gij*aT;
        C1 = 1;
        taut = (2*C1/dt)^2;
        taus = 1/sqrt((taut+tau1sqinv));
        TT(k) = 10*taus;
    end
    TT_ele(ele) = mean(TT);
end
end