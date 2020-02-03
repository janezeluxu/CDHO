function [u] = solveU(vertexData,BCb,IBC,BCval,meshData,solution_ien,...
    order_list,edge_usage,iper,tau,P)
%solve Kglobal*u = fGlobal matrix using gmres

global TotalDOF; 
global tol; 
[kGlobal, fGlobal] = globalKF(P,meshData,solution_ien,order_list,edge_usage,BCb,IBC,BCval,...
                              vertexData,iper,tau);
c = condest(kGlobal)                   
%u = gmres(kGlobal,fGlobal,50,tol,TotalDOF) ;  
display('linearsolveNow')
u = kGlobal\fGlobal;
end