function [u_n,ut_n] = time_step(vmax,vertexData,BCb,IBC,BCval,meshData,solution_ien,...
    order_ien,edge_usage,iper,OrderList,Area,IBC_vertex,uHBCnode)
global nstp;
global dt;
global Grid_size;
global node_Size;
global itau; 
[u_n,ut_n] = starting_u(IBC,BCval);
nstp
[JInverse, detJ,gij,x,y,A,hA] = getEleinfo(OrderList,meshData,solution_ien,order_ien,vertexData);
[patchArea,patchOrder] = getPatchArea(vertexData,order_ien,Area);
[c1,c2,PeHnode,tau] = initial_tau(vmax,vertexData,meshData,solution_ien,order_ien,Area,patchArea,patchOrder);
for istp = 1:nstp
    istp
    t = istp*dt;
    [BC_val,dBC_val] = eval_BC(IBC,BCval);
    [u_n,ut_n] = satisfy_boundary(IBC,BC_val,dBC_val,u_n,ut_n);

    [u_np1,ut_np1] = guess_flow(u_n,ut_n);

    [u_np1,ut_np1,tau_eleuse] = flow_solve(t,u_n,ut_n,u_np1,ut_np1,...
                     vertexData,BCb,IBC,BCval,meshData,solution_ien,...
                order_ien,edge_usage,iper,tau,OrderList,Area,IBC_vertex,...
                JInverse, detJ,gij);                      
    %% time update
    u_n = u_np1;
    ut_n = ut_np1;
    
    if mod(istp,31) == 0
%         display('writesolution');
%         % %% write the solution in vtk form
%         writeMesh(meshData,vertexData,'boundary')
%         solutionFile = strcat('./testvtk/cases/pulse/solution/p7dynamic',string(istp),'.txt');
%         writeSolution(u_np1,solutionFile)
%         writevtkfile(meshData,vertexData,u_n(1:node_Size),ut_n(1:node_Size),tau_eleuse);
%         output_filename = strcat('./testvtk/cases/pulse/dynamicp7/dynamic',string(istp),'.vtk')
%         twod_to_vtk_test (output_filename);
        
%         [taus,tauDC,taut_percent,tauadv_percent,taudiff_percent] = ...
%             getDC(t,u_n,ut_n,vertexData,JInverse,detJ,gij,meshData,solution_ien,order_ien,tau,OrderList);
%         writevtktau(taus,tauDC,taut_percent,tauadv_percent,taudiff_percent);
%         output_filename = strcat('./testvtk/cases/pulse/staticp1/static',string(istp),'tau.vtk');
%         twod_to_vtk_tau (output_filename);
    end
    
    if itau == 3            
        [tau] = solve_tau(t,u_n,ut_n,tau,...
                vertexData,meshData,solution_ien,order_ien, ...
                JInverse, detJ,gij,OrderList,IBC,Area,...
                patchArea,patchOrder,uHBCnode,istp,vmax);
    end
end

end

function [u_0,ut_0] = starting_u(IBC,BC_val)
global TotalDOF;
u_0 = zeros(TotalDOF,1);
ut_0 = zeros(TotalDOF,1);
for i = 1:length(IBC)
    u_0(IBC(i)) = BC_val(i);
    ut_0(IBC(i)) = 0;
end
end

function [u_np1,ut_np1] = guess_flow(u_n,ut_n)
global dt;
global gamma;
% ut_np1 = ut_n;
% %% use ODE1 to guess u_np1, p_np1
% u_np1 = u_n+dt*ut_n+dt*(ut_np1-ut_n);

%gamma = 1;
%% guess on ut_np1
ut_np1 = ((gamma-1)/gamma)*ut_n;
%% use ODE1 to guess u_np1, p_np1
u_np1 = u_n+dt*ut_n+gamma*dt*(ut_np1-ut_n);
end

function [u_n,ut_n] = satisfy_boundary(IBC,BC_val,dBC_val,u_n,ut_n)      
for i = 1:length(IBC)
    u_n(IBC(i)) = BC_val(i);
    ut_n(IBC(i)) = dBC_val(i);
end
end

function [BC_val,dBC_val] = eval_BC(IBC,BCval)
BC_val = BCval;
for i = 1:length(IBC)
    BC_val(i) = BC_val(i);
    dBC_val(i) = 0;
end
end

function [patchArea,patchOrder] = getPatchArea(vertexData,p_All,Area)
global node_Size;
patchArea = zeros(node_Size,1);
patchOrder = zeros(node_Size,1);
for node = 1:node_Size 
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
    mineleIEN = zeros(nElement,1);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        pAll = p_All(ele,:);
        patchArea(node) = patchArea(node)+Area(ele);
        mineleIEN(microEle) = min(pAll);
    end
    patchOrder(node) = min(mineleIEN);
end
end

function [c1node,c2node,PeHnode,tauele] = initial_tau(vmax,vertexData,...
                  meshData,solution_ien,p_ien,eleArea,patchArea,patchOrder)
 global Grid_size;
 global node_Size
 global kappa;
 global a;
c1node = zeros(node_Size,1);
c2node = zeros(node_Size,1);
PeHnode = zeros(node_Size,1);
tauele = ones(Grid_size,4)*5e-5;
 for node = 1:node_Size  
    c1node(node) = 1/vmax;
    
    AH = patchArea(node);
    const = 4/(sqrt(3)); 
    pPatch = patchOrder(node);
    
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
    
    velocity = [velocity1,velocity2];
                
    
    PeH = sqrt(const*AH)*norm(velocity)/(2*norm(kappa)*pPatch);
    PeHnode(node) = PeH;
    if PeH > 3
        c2node(node) = 1;
    else
        c2node(node) = 2;
    end
 end

 for ele = 1:Grid_size 
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    h = sqrt(eleArea(ele));
    %tauele = zeros(4,1);
    for i = 1:4
       node = IEN_mesh(i);
       tauele(ele,i) = c1node(node)*h^c2node(node);
    end
    %tau_ele(ele) = max(tauele);
 end
end