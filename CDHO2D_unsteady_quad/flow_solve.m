function [u_np1,ut_np1,tau_eleuse] = flow_solve(t,u_n,ut_n,u_np1,ut_np1,...
                     vertexData,BCb,IBC,BCval,meshData,solution_ien,...
                order_list,edge_usage,iper,tau,OrderList,Area,IBC_vertex,...
                JInverse, detJ,gij)
global Newton_tol;
global Newton_maxIter;
global Grid_size;
n = nIntergerPoints(OrderList,0);
qPoints = QuadGaussPoints(n);
nP = size(qPoints,1);
tau2 = zeros(Grid_size,nP);
tau1 = zeros(Grid_size,nP);

variable_start = u_n;
rhs_start = 1;
for iter = 1:Newton_maxIter
    iter
    [u,u_t] = iterPC(u_n,u_np1,ut_n,ut_np1);
    
    %[tau1,tau2,taut_percent,tauadv_percent,taudiff_percent] = ...
     %   getDC(t,u,u_t,JInverse,detJ,gij,meshData,solution_ien,order_list,tau,OrderList);
    [lhs, rhs,tau_eleuse] = globalKF(t,OrderList,meshData,solution_ien,order_list,edge_usage,BCb,IBC,BCval,...
        vertexData,iper,tau,u,u_t,tau1,tau2,JInverse, detJ,gij);
    
    [dut] = linear_solve(lhs,rhs);

    [u_np1,ut_np1] = update_flow(u_np1,ut_np1,dut);
    
    %% check convergence
    [ev_rhs,ev] = Newton_check(dut,rhs,variable_start,rhs_start)
    %rhs
    %convergence = [ev_rhs,ev]
    
%     if ev_rhs< Newton_tol(1) && ev<Newton_tol(2) 
%         break
%     end
    
end
end

function [u,u_t] =iterPC(u_n,u_np1,ut_n,ut_np1)
global af;
global am;
    %% solution at npaf and npam
    
    u = u_n+af*(u_np1-u_n);
    u_t = ut_n+am*(ut_np1-ut_n);

end

function [dvariable] = linear_solve(LHS,RHS)
dvariable = LHS\RHS;
end

function [ev_rhs,ev] = Newton_check (dvariable,rhs,variable_start,rhs_start,varargin)

if length(varargin) ==1
    vel_start = 1;
else
    vel_start = variable_start;
end
dvel = dvariable;

vres_start = rhs_start;

vres = rhs;


ev_rhs = norm(vres)/norm(vres_start);

ev = norm(dvel)/norm(vel_start);
end

function [u_np1,ut_np1] = update_flow(u_np1,ut_np1,dut)
global dt;
global gamma;

%% update ut, pt
ut_np1 = ut_np1+dut;

%% use ODE1 update u_np1, p_np1
u_np1 = u_np1+gamma*dt*dut;
end