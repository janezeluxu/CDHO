function [Error,H1] = errorEstimate(u,vertexData,meshData,solution_ien,p_ien,nInt,p)
        
global kappa; 
global casenumber;
global a;
global Grid_size;

%[ShapeFuncTable] = ShapeTable(nInt,p);
n = nInt;
qPoints = QuadGaussPoints(n);
[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);
%loop through all element
Error = 0;
H1_x = 0;
H1_y = 0;
for i = 1:Grid_size
    
    %[IENall,pAll] = elementIEN(i,IENstruct);
    [IEN_mesh,IENall,pAll] = elementIEN(i,meshData,solution_ien,p_ien);
    n = nInt;       
    qPoints = QuadGaussPoints(n);
    
    u_numerical = zeros(1,length(IENall));
    for j = 1:length(IENall)
        u_numerical(j) = u(IENall(j));
    end
    
    nP = size(qPoints,1);
    ux = zeros(nP,1);
    uy = zeros(nP,1);
    uquad = zeros(nP,1);
    
    xc = zeros(nP,1);
    yc = zeros(nP,1);
    
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
    
    for k = 1:nP
        xi = qPoints(k,1);
        eta = qPoints(k,2);
        [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
                    getjacobian(IEN_mesh,vertexData,xi,eta);

        uquad(k) = ShapeFunc(:,k)'*u_numerical';
        
        gradNGlobal = (divShapeFunc(:,:,k)*JInverse)';
        gradNx = gradNGlobal(1,:);
        gradNy = gradNGlobal(2,:);
        
        ux(k) = gradNx*u_numerical';
        uy(k) = gradNy*u_numerical';
        
        
        N1 = (1-xi)*(1-eta);
        N2 = xi*(1-eta);
        N3 = xi*eta;
        N4 = (1-xi)*eta;
        
        xc(k) = N1*x(1)+N2*x(2)+N3*x(3)+N4*x(4);
        yc(k) = N1*y(1)+N2*y(2)+N3*y(3)+N4*y(4);
    end
    
    
    [uA,ux_e,uy_e] = Exact(xc,yc,casenumber,a,kappa);
    %Error = (((uH - uA)).^2*qPoints(:,3)*detJ)+Error;
    Error = (((uquad' - uA')).^2*qPoints(:,3)*detJ)+Error;
    H1_x = (((ux' - ux_e')).^2*qPoints(:,3)*detJ)+H1_x;
    H1_y = (((uy' - uy_e')).^2*qPoints(:,3)*detJ)+H1_y;
    %if Error>1e-6
        %string = [ num2str(i),' Order is ', num2str(pAll), ' error is ',num2str(Error^0.5) ];
        %disp(string)
    %end
        
end
Error = Error^0.5;
H1 = (H1_x+H1_y)^0.5;
end





