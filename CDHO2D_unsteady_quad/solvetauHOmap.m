function [tauele,c1node,c2node,L,M] = solvetauHOmap(t,p,IBC,vertexData,meshData,...
    solution_ien,order_list,JInv, detJall,u,ut,patchArea,patchOrder,uHBCnode,eleArea,aref)

%global c2;
global kappa; 
global a;
global node_Size;
global Grid_size;
global casenumber;
tauele = zeros(Grid_size,4);
c1node = zeros(node_Size,1);
c2node = zeros(node_Size,1);
L = zeros(node_Size,2);
M = zeros(node_Size,2);
%Hhfactor = getHhfactor(vertexData,node_Size);
%nInt = 0;
%n = 2*p;
n = nIntergerPoints(p,0);
qPoints = QuadGaussPoints(n); 
[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);
xpower = zeros(p+2,1);
ypower = zeros(p+2,1);
for node = 1:node_Size 
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
    %PeHDis(node) = PeH;
    if PeH > 3
        c2 = 1;
        c2node(node) = 1;
    else
        c2 = 2;
        c2node(node) = 2;
    end
        
    SurroundEle = vertexData{node,3};
    nElement = length(SurroundEle);    
    %% get LS intergral        
    %n = 2*p;
    %qPoints = QuadGaussPoints(n);  
    xLS = zeros(size(qPoints,1),nElement);
    yLS = zeros(size(qPoints,1),nElement);
    uLS = zeros(size(qPoints,1),nElement);
    utLS = zeros(size(qPoints,1),nElement);
    w = qPoints(:,3);
    
    J = zeros(1,nElement);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
        
        %% get LS        
        %[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,pAll(1));
        uA = zeros(1,length(IENall));
        uAt = zeros(1,length(IENall));
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
            uAt(j) = ut(IENall(j));
        end
        uh = uA*ShapeFunc;
        uht = uAt*ShapeFunc;
        
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

        detJele = detJall(ele,:);
        for k = 1:nP
            xi = qPoints(k,1);
            eta = qPoints(k,2);
            detJ = detJele(k);
                
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
        utLS(:,microEle) = uht';
        J(microEle) = detJ;
    end
    pH = max(min(pAll),1);
    xmin = min(min(xLS));
    xmax=max(max(xLS));
    ymin=min(min(yLS));
    ymax=max(max(yLS));
    LS_u = getLSInt(xLS,yLS,xmin,xmax,ymin,ymax,uLS,pH,uHBCnode(node,:),vertexData,nElement,w,J);
    LS_ut = getLSInt(xLS,yLS,xmin,xmax,ymin,ymax,utLS,pH,uHBCnode(node,:),vertexData,nElement,w,J);
              
    BB=0;CC=0;DD=0;
    AA=0;EE=0;FF=0;
    %Hh = Hhfactor(node);
    
    H = sqrt(AH);
    for microEle = 1:nElement
        ele = SurroundEle(microEle);
        h = sqrt(eleArea(ele));
        [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);   
        %% get uh intergrals
        %n = nIntergerPoints(max(pAll),0);
        %qPoints = QuadGaussPoints(n);
        
         uA = zeros(1,length(IENall));
         uAt = zeros(1,length(IENall));
         
        %[ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,pAll(1));
 
        for j = 1:length(IENall)
            uA(j) = u(IENall(j));
            uAt(j) = ut(IENall(j));
        end
        
        vIDs = IEN_mesh;
        x1 = vertexData{vIDs(1),2}(1);
        y1 = vertexData{vIDs(1),2}(2);
        x2 = vertexData{vIDs(2),2}(1);
        y2 = vertexData{vIDs(2),2}(2);
        x3 = vertexData{vIDs(3),2}(1);
        y3 = vertexData{vIDs(3),2}(2);
        x4 = vertexData{vIDs(4),2}(1);
        y4 = vertexData{vIDs(4),2}(2);
        
        xCord = [x1,x2,x3,x4];
        yCord = [y1,y2,y3,y4];
            
        JInvele = JInv(ele,:);
        detJele = detJall(ele,:);
        %gijele = gijall(ele,:);
        count = 1;
        for i = 1:length(qPoints) 
            
            xi = qPoints(i,1);
            eta = qPoints(i,2);
            
            JInverse = reshape(JInvele(count:count+3),2,2);
            detJ = detJele(k);
            %gij = reshape(gijele(count:count+3),2,2);
            count = count+4;
            
            Jw = detJ*qPoints(i,3);
            N1 = (1-xi)*(1-eta);
            N2 = xi*(1-eta);
            N3 = xi*eta;
            N4 = (1-xi)*eta;
            
            x = N1*xCord(1)+N2*xCord(2)+N3*xCord(3)+N4*xCord(4);
            y = N1*yCord(1)+N2*yCord(2)+N3*yCord(3)+N4*yCord(4);
           
            f = MMS(t,x,y,casenumber);
            [uH,duH]=calcuH(x,y,LS_u,pH,xmin,xmax,ymin,ymax,node,xpower,ypower);
            [uHt,duHt]=calcuH(x,y,LS_ut,pH,xmin,xmax,ymin,ymax,node,xpower,ypower);
            
            WHx = [duH(1),0;0,duH(2)];
            WH = uH;
            %WHx = [1,0;0,1];
            
            if isa(a{1}, 'function_handle')==1
                velocity1 = a{1}(x,y);
            else
                velocity1 = a{1};
            end
            if isa(a{2}, 'function_handle')==1
                velocity2 = a{2}(x,y);
            else
                velocity2 = a{2};
            end
            
            velocity = [velocity1,velocity2];
    
            uh = uA*ShapeFunc(:,i);
            duh = uA*divShapeFunc(:,:,i)*JInverse;
            uht = uAt*ShapeFunc(:,i);
            %% WH1 = [WHX,0];WH2 = [0,WHY];
%             BB = BB+WHx*H^(c2)*velocity'*(uHt+velocity*duH-f)*Jw;
%             CC = CC+WH*[velocity(1)*duH(1);velocity(2)*duH(2)]*Jw+WH*uHt*Jw;%-WH*f*Jw;
%             DD = DD+WHx*kappa*duH*Jw;
%             
%             AA = AA+WHx*h^(c2)*velocity'*(uht+velocity*duh'-f)*Jw;
%             EE = EE+WH*[velocity(1)*duh(1);velocity(2)*duh(2)]*Jw+WH*uht*Jw;%-WH*f*Jw;
%             FF = FF + WHx*kappa*duh'*Jw;
            
              %% WH choice2
            BB = BB+H^(c2)*(WHx(1)*velocity(1)+WHx(2)*velocity(2))*(uHt+velocity*duH-f)*Jw;
            CC = CC+WH*(velocity(1)*duH(1)+velocity(2)*duH(2))*Jw+WH*uHt*Jw;%-WH*f*Jw;
            
            AA = AA+h^(c2)*(WHx(1)*velocity(1)+WHx(2)*velocity(2))*(uht+velocity*duh'-f)*Jw;
            EE = EE+WH*(velocity(1)*duh(1)+velocity(2)*duh(2))*Jw+WH*uht*Jw;%-WH*f*Jw;
            
%               %% WH choice3
%             BB = BB+WHx*H^(c2)*velocity'*(uHt+velocity*duH-f)*Jw;
%             CC = CC-uH*[velocity(1)*WHx(1);velocity(2)*WHx(2)]*Jw;%+WH*uHt*Jw;%-WH*f*Jw;
%             
%             AA = AA+WHx*h^(c2)*velocity'*(uht+velocity*duh'-f)*Jw;
%             EE = EE-uh*[velocity(1)*WHx(1);velocity(2)*WHx(2)]*Jw;%+WH*uht*Jw;%-WH*f*Jw;
 
%               %% WH choice4
%             BB = BB+H^(c2)*(WHx(1)*velocity(1)+WHx(2)*velocity(2))*(uHt+velocity*duH-f)*Jw;
%             CC = CC-uH*(velocity(1)*WHx(1)+velocity(2)*WHx(2))*Jw+WH*uHt*Jw;%-WH*f*Jw;
%             
%             AA = AA+h^(c2)*(WHx(1)*velocity(1)+WHx(2)*velocity(2))*(uht+velocity*duh'-f)*Jw;
%             EE = EE-uh*(velocity(1)*WHx(1)+velocity(2)*WHx(2))*Jw+WH*uht*Jw;%-WH*f*Jw;
        end
    end
    M(node,:) = AA-BB;
    L(node,:) = (EE-CC);
    
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

%% pathtube average
Lave = zeros(node_Size,1);
Mave = zeros(node_Size,1);

for node = 1:node_Size 
    SurroundNode = vertexData{node,4};
    nNode = length(SurroundNode);
    Lg = 0; Mg = 0;
    for nN = 1:nNode
        micronode = SurroundNode(nN);
        Lg = Lg+abs(L(micronode));
        Mg = Mg+abs(M(micronode));
    end
    Lave(node) = Lg/nNode;
    Mave(node) = Mg/nNode;
end

% for node = 1:node_Size 
%     if norm(M(node,:))<1e-10
%        c1node(node) = 0;
%     else
%         MM = M(node,:);
%         LL = L(node,:);
%         %taunodenode=(MM')\(LL');
%         c1node(node) = (MM')\(LL');
%     end
% end

for node = 1:node_Size 
    if norm(M(node,:))<1e-10
       c1node(node) = 0;
    else
        MM = Mave(node);
        LL = Lave(node);
        %taunodenode=(MM')\(LL');
        c1node(node) = (MM)\(LL);
    end
end

% for node = 1:node_Size   
%     if c1node(node) < 0
%         %c1node(node)=0;
%         c1node(node)=abs(c1node(node));
%     end
%     
% %     if c1node(node) > 10
% %         %c1node(node)=0;
% %         c1node(node) = 10;
% %     end
% end

c1node_S = std(c1node)
c1node_mean = mean(c1node)

for node = 1:node_Size   
    if c2node(node)>1
        if c1node(node) > 1/norm(kappa)/p%c1node_mean+c1node_S
        %c1node(node)=0;
        c1node(node) = 1/norm(kappa)/p;%c1node_mean+c1node_S;
        end
    else
    if c1node(node) > 1/aref/p%c1node_mean+c1node_S
        %c1node(node)=0;
        c1node(node) = 1/aref/p;%c1node_mean+c1node_S;
    end
    end
end

for i = 1:length(IBC)
    c1node(IBC(i)) = 0;
end

for ele = 1:Grid_size 
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,order_list);
    h = sqrt(eleArea(ele));
    for i = 1:4
        node = IEN_mesh(i);
        tauele(ele,i) = c1node(node)*h^c2node(node);
    end
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
%lhs = zeros(sizeA,sizeA);
%rhs = zeros(sizeA,1);

% r=0;
% for k1 = 0:p
%     for m1 = 0:k1
%         r=r+1;
%         c=0;
%         for k2 = 0:p
%             for m2 = 0:k2
%                 c=c+1;
%                     sumxyn = 0;
%                     for l = 1:size(x,1)
%                         for ele = 1:nEle
%                             sumxyn = sumxyn+x(l,ele)^(k1-m1)*y(l,ele)^(m1)*...
%                             x(l,ele)^(k2-m2)*y(l,ele)^(m2)*J(ele)*w(l);
%                         end
%                     end
%                 lhs(r,c) = sumxyn;
%             end
%         end
%     end
% end

a = zeros(sizeA,sizeA);
b = zeros(sizeA,sizeA);
r=0;
for k1 = 0:p
    for m1 = 0:k1
        r=r+1;
        c=0;
        for k2 = 0:p
            for m2 = 0:k2
                c=c+1;               
                a(r,c) = k1-m1+k2-m2;
                b(r,c) = m1+m2;
            end
        end
    end
end

x1size = size(x,1);
Jw = zeros(x1size,nEle);
for l = 1:size(x,1)
    for ele = 1:nEle
        Jw(l,ele) = J(ele)*w(l);
    end
end

%lhs1 = zeros(sizeA,sizeA);
% for i = 1:sizeA
%     for j = i:sizeA
%         lhs1(i,j) = sum(sum(x.^(a(i,j)).*y.^(b(i,j)).*Jw));
%     end
% end
% lhs1 = lhs1+lhs1'-diag(diag(lhs1)); 
%lhstest-lhs1

maxa = max(max(a));
maxb = max(max(b));
xpower = zeros(x1size,nEle,maxa);
ypower = zeros(x1size,nEle,maxb);

for i = 0:maxa
    xpower(:,:,i+1) = x.^i;
end

for i = 0:maxb
    ypower(:,:,i+1) = y.^i;
end

lhs1 = zeros(sizeA,sizeA);
for i = 1:sizeA
    for j = i:sizeA
        xp = xpower(:,:,a(i,j)+1);
        yp = ypower(:,:,b(i,j)+1);
        lhs1(i,j) = sum(sum(xp.*yp.*Jw));
    end
end
lhs1 = lhs1+lhs1'-diag(diag(lhs1)); 
%lhs1test-lhs1

rhs = zeros(sizeA,1);
%rhs = zeros(sizeA,1);
c=0;
for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        xp = xpower(:,:,(k1-m1)+1);
        yp = ypower(:,:,m1+1);
        rhs(c) = sum(sum(u.*xp.*yp.*Jw));
%        rhs(c) = sum(sum(u.*x.^(k1-m1).*y.^(m1).*Jw));
%         sumeleu = 0;       
%         for l = 1:size(x,1)
%             for ele = 1:nEle
%                 sumeleu = sumeleu+u(l,ele)*x(l,ele)^(k1-m1)*y(l,ele)^(m1)*J(ele)*w(l);
%             end
%         end
%        rhs(c) = sumeleu;
    end
end
%rhs-rhs

LS = lhs1\rhs;

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

function [uH,duH]=calcuH(x,y,LS,p,xmin,xmax,ymin,ymax,node,xpower,ypower)
%% mapping to a box
[xi,eta]= mapxy(x,y,xmin,xmax,ymin,ymax);
[uH,duHxi,duHeta]=getuH(xi,eta,LS,p,xpower,ypower);
duHx=duHxi/(xmax-xmin);
duHy = duHeta/(ymax-ymin);
duH=[duHx;duHy];
%[uH,duHx,duHy]=getuH(x,y,LS,p);
%duH=[duHx;duHy];
end

function [uHtest,duHxtest,duHytest]=getuH(x,y,LS,p,xpower,ypower)
%% uH = a0+a1x+a2y+a3x^2+a4xy+a5y^2+...
% duHx = 0;
% duHy = 0;
% c=0;
% uH = 0;
% 
% for k1 = 0:p
%     for m1 = 0:k1
%         c=c+1;
%         uH = uH+LS(c)*x^(k1-m1)*y^(m1);
%         if (k1-m1) ==0
%             addduHx=0;
%         else
%             addduHx = (k1-m1)*LS(c)*x^(k1-m1-1)*y^(m1);
%         end
%         %duHx = duHx+(k1-m1)*LS(c).*x.^(k1-m1-1).*y.^(m1);
%         duHx = duHx+addduHx;
%         if (m1) == 0
%             addduHy=0;
%         else
%             addduHy = (m1)*LS(c)*x^(k1-m1)*y^(m1-1);
%         end
%         duHy = duHy+addduHy;
%         %duHy = duHy+(m1)*LS(c).*x.^(k1-m1).*y.^(m1-1);
%     end
% end


duHxtest = 0;
duHytest = 0;
c=0;
uHtest = 0;

%xpower = zeros(p+2,1);
%ypower = zeros(p+2,1);
for i = 0:p
    xpower(i+2) = x^i;
    ypower(i+2) = y^i;
end

%addduHx=0;
%addduHy=0;
for k1 = 0:p
    for m1 = 0:k1
        c=c+1;
        uHtest = uHtest+LS(c)*xpower(k1-m1+2)*ypower(m1+2);
        %addduHx = (k1-m1)*LS(c)*xpower(k1-m1+1)*ypower(m1+2);
        duHxtest = duHxtest+(k1-m1)*LS(c)*xpower(k1-m1+1)*ypower(m1+2);
        %addduHy = (m1)*LS(c)*xpower(k1-m1+2)*ypower(m1+1);
        duHytest = duHytest+(m1)*LS(c)*xpower(k1-m1+2)*ypower(m1+1);
    end
end
% uHtest-uH
% duHxtest-duHx
% duHytest-duHy

% duHxtest = 0;
% duHytest = 0;
% c=0;
% uHtest = 0;
% 
% Index1 = zeros(length(LS),1);
% Index2 = zeros(length(LS),1);
% for k1 = 0:p
%     for m1 = 0:k1
%         c=c+1;
%         Index1(c) = (k1-m1);
%         Index2(c) = m1;
%     end
% end
% 
% addduHx = 0;
% addduHy=0;
% for c = 1:length(LS)
%         uHtest = uHtest+LS(c)*x^(Index1(c))*y^(Index2(c));
%         %if (k1-m1) ==0
%         %    addduHx=0;
%         if (k1-m1) >0
%             addduHx = Index1(c)*LS(c)*x^(Index1(c)-1)*y^(Index2(c));
%         end
%         %duHx = duHx+(k1-m1)*LS(c).*x.^(k1-m1-1).*y.^(m1);
%         duHxtest = duHxtest+addduHx;
%         %if (m1) == 0
%         %    addduHy=0;
%         if (m1) > 0
%             addduHy = Index2(c)*LS(c)*x^(Index1(c))*y^(Index2(c)-1);
%         end
%         duHytest = duHytest+addduHy;
%         %duHy = duHy+(m1)*LS(c).*x.^(k1-m1).*y.^(m1-1);
% end
end
