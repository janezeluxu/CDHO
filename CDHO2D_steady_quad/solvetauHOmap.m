function [taunode,L,M] = solvetauHOmap(p,vertexData,patchArea,patchOrder,IBC,...
    uHBCnode,meshData,solution_ien,order_list,u)
%global c2;
global kappa; 
global a;
global node_Size;
global Grid_size;
n = nIntergerPoints(p,0);
qPoints = QuadGaussPoints(n);
[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);

taunode = zeros(node_Size,1);
c2Dis = zeros(node_Size,1);
L = zeros(node_Size,2);
M = zeros(node_Size,2);
Hhfactor = getHhfactor(vertexData,node_Size);
nInt = 0;
for node = 1:node_Size 
    AH = patchArea(node);
    const = 4/(sqrt(3)); 
    pPatch = patchOrder(node);
    PeH = sqrt(const*AH)*norm(a)/(2*norm(kappa)*pPatch);
    %PeHDis(node) = PeH;
    if PeH > 3
        c2 = 1;
        c2Dis(node) = 1;
    else
        c2 = 2;
        c2Dis(node) = 2;
    end
        
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);
%     xLS = [];yLS=[];uLS=[];
%     for microEle = 1:nElement
%         ele = SurroundEle(microEle);
%         [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
%     
%         uA = zeros(1,length(IENall));
% 
%         for j = 1:length(IENall)
%             uA(j) = u(IENall(j));
%         end
%         %uh = uA*ShapeFunc
%         
%         vIDs = IEN_mesh;
%         x1 = vertexData{vIDs(1),2}(1);
%         y1 = vertexData{vIDs(1),2}(2);
%         x2 = vertexData{vIDs(2),2}(1);
%         y2 = vertexData{vIDs(2),2}(2);
%         x3 = vertexData{vIDs(3),2}(1);
%         y3 = vertexData{vIDs(3),2}(2);
%         x4 = vertexData{vIDs(4),2}(1);
%         y4 = vertexData{vIDs(4),2}(2);
%         
%         xcord = [x1,x2,x3,x4];
%         ycord = [y1,y2,y3,y4];
%         
%         dx = 0:1/p:1;
%         count = 1;
%         for a1 = 1:length(dx)
%             for a2 = 1:length(dx)
%                 eval_points(count,1) = dx(a1);
%                 eval_points(count,2) = dx(a2);
%                 count= count+1;
%             end
%         end
%        nP = length(eval_points);
%        x = zeros(1,nP);
%         y = zeros(1,nP);
%         uh = zeros(1,nP);
%         
%        [ ShapeFunc, ~ ] = ShapeFunc_quad(eval_points,p);
%         for k = 1:nP
%             xi = eval_points(k,1);
%             eta = eval_points(k,2);
%             N1 = (1-xi)*(1-eta);
%             N2 = xi*(1-eta);
%             N3 = xi*eta;
%             N4 = (1-xi)*eta;
%             x(k) = N1*xcord(1)+N2*xcord(2)+N3*xcord(3)+N4*xcord(4);
%             y(k) = N1*ycord(1)+N2*ycord(2)+N3*ycord(3)+N4*ycord(4);
%             
%             uh(k) = uA*ShapeFunc(:,k);
%         end
%         
%         xLS = [xLS,x];
%         yLS = [yLS,y];
%         uLS = [uLS,uh];
% 
%     end
% 
%     pH = max(min(pAll)-1,1);
%     xmin = min(xLS); xmax=max(xLS); ymin=min(yLS);ymax=max(yLS);
%     %[xi,eta]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax);
%     LS=getLS(xLS,yLS,xmin,xmax,ymin,ymax,uLS,pH,uHBCnode(node,:),vertexData);
%     
    %% get LS intergral        
    n = 2*p;
    qPoints = QuadGaussPoints(n);  
    xLS = zeros(size(qPoints,1),nElement);
    yLS = zeros(size(qPoints,1),nElement);
    uLS = zeros(size(qPoints,1),nElement);
    w = qPoints(:,3);
    
    J = zeros(1,nElement);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
        
        %% get LS        
        [ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,pAll(1));
        uA = zeros(1,length(IENall));
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
        end
        uh = uA*ShapeFunc;
        
        nP = length(qPoints);
        x = zeros(1,nP);
        y = zeros(1,nP);
        
        vIDs = IEN_mesh;
        x1 = vertexData{vIDs(1),2}(1);
        y1 = vertexData{vIDs(1),2}(2);
        x2 = vertexData{vIDs(2),2}(1);
        y2 = vertexData{vIDs(2),2}(2);
        x3 = vertexData{vIDs(3),2}(1);
        y3 = vertexData{vIDs(3),2}(2);
        x4 = vertexData{vIDs(4),2}(1);
        y4 = vertexData{vIDs(4),2}(2);
        
        xcord = [x1,x2,x3,x4];
        ycord = [y1,y2,y3,y4];

        for k = 1:nP
            xi = qPoints(k,1);
            eta = qPoints(k,2);
            [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                getjacobian(IEN_mesh,vertexData,xi,eta);
            N1 = (1-xi)*(1-eta);
            N2 = xi*(1-eta);
            N3 = xi*eta;
            N4 = (1-xi)*eta;
            x(k) = N1*xcord(1)+N2*xcord(2)+N3*xcord(3)+N4*xcord(4);
            y(k) = N1*ycord(1)+N2*ycord(2)+N3*ycord(3)+N4*ycord(4);
            
        end
        
        xLS(:,microEle) = x';
        yLS(:,microEle) = y';
        uLS(:,microEle) = uh';
        J(microEle) = detJ;
    end
    pH = max(min(pAll),1);
    xmin = min(min(xLS));
    xmax=max(max(xLS));
    ymin=min(min(yLS));
    ymax=max(max(yLS));
    LS = getLSInt(xLS,yLS,xmin,xmax,ymin,ymax,uLS,pH,uHBCnode(node,:),vertexData,nElement,w,J);

    
    BB=0;CC=0;DD=0;
    AA=0;EE=0;FF=0;
    Hh = Hhfactor(node);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);   
        
        %% get uH intergrals
        n = nIntergerPoints(pH,nInt);
        qPoints = QuadGaussPoints(n);
        
%         for i = 1:length(qPoints)
%             
%             xi = qPoints(i,1);
%             eta = qPoints(i,2);
%             [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
%                 getjacobian(IEN_mesh,vertexData,xi,eta);
%             Jw = detJ*qPoints(i,3);
%             N1 = (1-xi)*(1-eta);
%             N2 = xi*(1-eta);
%             N3 = xi*eta;
%             N4 = (1-xi)*eta;
%             
%             x = N1*xCord(1)+N2*xCord(2)+N3*xCord(3)+N4*xCord(4);
%             y = N1*yCord(1)+N2*yCord(2)+N3*yCord(3)+N4*yCord(4);
%            
%             [uH,duH]=calcuH(x,y,LS,pH,xmin,xmax,ymin,ymax,node);
%             BB = BB+Hh^(c2)*(a'*a)*duH*Jw;
%             CC = CC+a'*uH*Jw;
%             DD = DD+kappa*duH*Jw;
%             
%         end
        %% get uh intergrals
        n = nIntergerPoints(max(pAll),0);
        qPoints = QuadGaussPoints(n);
        
         uA = zeros(1,length(IENall));
         
         [ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,pAll(1));
 
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
        end
        for i = 1:length(qPoints) 
            
            xi = qPoints(i,1);
            eta = qPoints(i,2);
            [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                getjacobian(IEN_mesh,vertexData,xi,eta);
            Jw = detJ*qPoints(i,3);
            N1 = (1-xi)*(1-eta);
            N2 = xi*(1-eta);
            N3 = xi*eta;
            N4 = (1-xi)*eta;
            
            x = N1*xCord(1)+N2*xCord(2)+N3*xCord(3)+N4*xCord(4);
            y = N1*yCord(1)+N2*yCord(2)+N3*yCord(3)+N4*yCord(4);
           
            [uH,duH]=calcuH(x,y,LS,pH,xmin,xmax,ymin,ymax,node);
            WHx = [duH(1),0;0,duH(2)];
            %WHx = [1,0;0,1];
            
            BB = BB+WHx*Hh^(c2)*(a'*a)*duH*Jw;
            CC = CC+WHx*a'*uH*Jw;
            DD = DD+WHx*kappa*duH*Jw;
            
%             xi = qPoints(i,1);
%             eta = qPoints(i,2);
%             [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
%                 getjacobian(IEN_mesh,vertexData,xi,eta);
%             Jw = detJ*qPoints(:,3);
            uh = uA*ShapeFunc(:,i);
            duh = uA*divShapeFunc(:,:,i)*JInverse;
        
            AA = AA+WHx*(a'*a)*duh'*Jw;
            EE = EE+WHx*a'*uh*Jw;
            FF = FF + WHx*kappa*duh'*Jw;
        end
    end
    M(node,:) = AA-BB;
    L(node,:) = EE-CC;
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if norm(M)<1e-20
%         taunode(node) = 0;
%     else
%         taunode(node) = M\L;
%     end
%     
%     %taunode(node) = M\L;
%      
%      if taunode(node) < 0
%          taunode(node)=0;
%          %taunode(node)=abs(taunode(node));
%      end
end

%aref = norm(a);
for node = 1:node_Size 
    if norm(M(node,:))<1e-10
        taunode(node) = 0;
    else
        MM = M(node,:);
        LL = L(node,:);
        %taunodenode=(MM')\(LL');
        taunode(node) = (MM')\(LL');
    end
    
    if taunode(node) < 0
        taunode(node)=0;
        %taunode(node)=abs(taunode(node));
    end
    
%     if c2Dis(node)>1
%         countc2 = countc2+1;
%         if taunode(node) > 10/kappa/p/p
%             taunode(node) = 10/kappa/p/p;
%             %count1  = count1+1;
%         end
%         
%         if taunode(node) < 0.001/kappa/p/p
%             taunode(node) = 0.001/kappa/p/p;
%             %count2  = count2+1;
%         end
%     else
%         if taunode(node) > 10/aref/p
%             taunode(node) = 10/aref/p;
%             %count3  = count3+1;
%         end
%         
%         if taunode(node) < 0.001/aref/p
%            %taunode(node) 
%            taunode(node) = 0.001/aref/p;
%            %count4  = count4+1;
%         end
%     end
end

% %% pathtube average
% LM = zeros(node_Size,1);
% MM = zeros(node_Size,1);
% for i = 1:node_Size
%     LM(i) = (L(i,1)*M(i,1)+L(i,2)*M(i,2));
%     MM(i) = (M(i,1)*M(i,1)+M(i,2)*M(i,2));
%     %testLM(i) = MM(i)\LM(i);
% end
% 
% for node = 1:node_Size 
%     SurroundNode = vertexData{node,4};
%     nNode = length(SurroundNode);
%     Lg = [0,0]; Mg = [0,0];
%     for nN = 1:nNode
%         micronode = SurroundNode(nN);
%         Lg = Lg+L(micronode,:);
%         Mg = Mg+M(micronode,:);
%     end
%     
%     if norm(Mg)<1e-16
%         taunode(node) = 0;
%     else
%         taunode(node) = (Mg')\(Lg');
%     end
%     
%     if taunode(node) < 0
%         %taunode(node)=0;
%         taunode(node)=abs(taunode(node));
%     end
% end

% for i = 1:length(IBC)
%     j = IBC(i);
%     taunode(j) = 0;
% end
end

function [xi,eta]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax)
xi = (xLS-xmin)./(xmax-xmin);
eta = (yLS-ymin)./(ymax-ymin);
end

function LS = getLSInt(xLS,yLS,xmin,xmax,ymin,ymax,u,p,uHBCnode,vertexData,nEle,w,J)
[x,y]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax);

sizeA = (p+1)*(p+2)/2;
lhs = zeros(sizeA,sizeA);
rhs = zeros(sizeA,1);

r=0;
for k1 = 0:p
    for m1 = 0:k1
        r=r+1;
        c=0;
        for k2 = 0:p
            for m2 = 0:k2
                c=c+1;
                    sumxyn = 0;
                    for l = 1:size(x,1)
                        for ele = 1:nEle
                            sumxyn = sumxyn+x(l,ele)^(k1-m1)*y(l,ele)^(m1)*...
                            x(l,ele)^(k2-m2)*y(l,ele)^(m2)*J(ele)*w(l);
                        end
                    end
                lhs(r,c) = sumxyn;
            end
        end
    end
end

c=0;
for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        sumeleu = 0;
        for l = 1:size(x,1)
            for ele = 1:nEle
                sumeleu = sumeleu+u(l,ele)*x(l,ele)^(k1-m1)*y(l,ele)^(m1)*J(ele)*w(l);
            end
        end
        rhs(c) = sumeleu;
    end
end

LS = lhs\rhs;

flag = uHBCnode(1);
bcNode = uHBCnode(2);
bcValue = uHBCnode(3);
if flag
    %% exact fit bcNode
    %a = bcValue-LS(2)*xx(end)-LS(3)*xx(end)*xx(end);
    atest = 0;
    c=0;
    xbcnode = vertexData{bcNode,2}(1);
    ybcnode = vertexData{bcNode,2}(2);
    [xbc,ybc]=mapxy(xbcnode,ybcnode,xmin,xmax,ymin,ymax);
    
    %% need mapping
    for k1 = 0:p
        for m1 = 0:k1
            c=c+1;
            atest = atest-LS(c).*xbc^(k1-m1).*ybc^(m1);
        end
    end
    
    a = atest+bcValue+LS(1);
    LS(1) = a;    
end

%LSIntflag = LS
end

function LS = getLS(xLS,yLS,xmin,xmax,ymin,ymax,u,p,uHBCnode,vertexData)
[x,y]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax);
%x = xLS;
%y = yLS;
sizeA = (p+1)*(p+2)/2;
lhs = zeros(sizeA,sizeA);
rhs = zeros(sizeA,1);

xyn = zeros(1,length(x));
ufxn = zeros(1,length(x));

r=0;
for k1 = 0:p
    for m1 = 0:k1
        r=r+1;
        c=0;
        for k2 = 0:p
            for m2 = 0:k2
                c=c+1;
                    for l = 1:length(x)             
                        xyn(l) = x(l)^(k1-m1)*y(l)^(m1)*x(l)^(k2-m2)*y(l)^(m2);
                    end
                    lhs(r,c) = sum(xyn);
            end
        end
    end
end
c=0;
for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        for l = 1:length(x)
            ufxn(l) = u(l)*x(l)^(k1-m1)*y(l)^(m1);
        end
        rhs(c) = sum(ufxn);
    end
end
LS = lhs\rhs;

flag = uHBCnode(1);
bcNode = uHBCnode(2);
bcValue = uHBCnode(3);
if flag 
   %% exact fit bcNode
   %a = bcValue-LS(2)*xx(end)-LS(3)*xx(end)*xx(end);
   atest = 0;
   c=0;
   xbcnode = vertexData{bcNode,2}(1);
   ybcnode = vertexData{bcNode,2}(2);
   [xbc,ybc]=mapxy(xbcnode,ybcnode,xmin,xmax,ymin,ymax);
   
   %% need mapping   
   for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        atest = atest-LS(c).*xbc^(k1-m1).*ybc^(m1);
    end
   end
   
   a = atest+bcValue+LS(1);
   LS(1) = a;
   
end

end

function [uH,duH]=calcuH(x,y,LS,p,xmin,xmax,ymin,ymax,node)
%% mapping to a box
[xi,eta]= mapxy(x,y,xmin,xmax,ymin,ymax);
[uH,duHxi,duHeta]=getuH(xi,eta,LS,p);
duHx=duHxi/(xmax-xmin);
duHy = duHeta/(ymax-ymin);
duH=[duHx;duHy];
%[uH,duHx,duHy]=getuH(x,y,LS,p);
%duH=[duHx;duHy];
end

function [uH,duHx,duHy]=getuH(x,y,LS,p)
%% uH = a0+a1x+a2y+a3x^2+a4xy+a5y^2+...
duHx = 0;
duHy = 0;
c=0;
uH = 0;
for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        uH = uH+LS(c).*x.^(k1-m1).*y.^(m1);
        if (k1-m1) ==0
            addduHx=0;
        else
            addduHx = (k1-m1)*LS(c).*x.^(k1-m1-1).*y.^(m1);
        end
        %duHx = duHx+(k1-m1)*LS(c).*x.^(k1-m1-1).*y.^(m1);
        duHx = duHx+addduHx;
        if (m1) == 0
            addduHy=0;
        else
            addduHy = (m1)*LS(c).*x.^(k1-m1).*y.^(m1-1);
        end
        duHy = duHy+addduHy;
        %duHy = duHy+(m1)*LS(c).*x.^(k1-m1).*y.^(m1-1);
    end
end
end

function Hh = getHhfactor(vertexData,TotalNode)
Hh = zeros(TotalNode,1);
for node = 1:TotalNode 
    SurroundEle = vertexData{node,3};
    Hh(node) = sqrt(length(SurroundEle));
end
Hh = 2*ones(TotalNode,1);
end