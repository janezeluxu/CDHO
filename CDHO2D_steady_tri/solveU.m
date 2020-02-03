function [u] = solveU(vertexData,BCb,IBC,BCval,meshData,solution_ien,...
    order_list,edge_usage,iper,tau,varargin)
%solve Kglobal*u = fGlobal matrix using gmres

global TotalDOF; 
global tol; 


if length(varargin) ==1
    OrderList = varargin{1};
    P = OrderList;
elseif length(varargin) ==3
    PV = varargin{1};
    PH = varargin{2};
    PT = varargin{3};
    P = [PV,PH,PT];
    P = unique(P);
end

[kGlobal, fGlobal] = globalKF(P,meshData,solution_ien,order_list,edge_usage,BCb,IBC,BCval,...
                              vertexData,iper,tau);
c = condest(kGlobal)                   
%u = gmres(kGlobal,fGlobal,50,tol,TotalDOF) ;  
display('linearsolveNow')
u = kGlobal\fGlobal;

% slavenodes = [5   122   123     6   124   125     7   126   127     8   ,...
%     128   129     9   130   131    10,...
% 132   133    11   134   135    12   136   137    13   138   139     3];
% 
% masternodes = [31   181   180    30   179   178    28   177   176    29   175 ,...
%     174    23   173   172    27, 171   170    26   169   168    24   167  ,...
%     166    25   165   164     4];
% for i = 1:length(slavenodes)
%     slave = u(slavenodes(i))
%     master = u(masternodes(i))
% end
end