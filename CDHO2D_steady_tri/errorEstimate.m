function [Error,H1] = errorEstimate(u,vertexData,meshData,solution_ien,p_ien,nInt,varargin)
        
global kappa; 
global force;
global a;
global Grid_size;
if length(varargin) ==1
    OrderList = varargin{1};
    p = OrderList;
elseif length(varargin) ==3
    PV = varargin{1};
    PH = varargin{2};
    PT = varargin{3};
    p = [PV,PH,PT];
end

[ShapeFuncTable] = ShapeTable(nInt,p);

%loop through all element
Error = 0;
H1_x = 0;
H1_y = 0;
for i = 1:Grid_size
    
    %[IENall,pAll] = elementIEN(i,IENstruct);
    [IEN_mesh,IENall,pAll] = elementIEN(i,meshData,solution_ien,p_ien);
    uE = zeros(1,length(IENall));
    n = nInt;       
    qPoints = TriGaussPoints(n);
        
    if range(pAll) == 0
        %It is a uniform element 
        %ShapeFunc = ShapeFuncTable{pAll(1)};
        simplexsf = SimplexShapeFunc(qPoints,pAll(1));
        [ShapeFunc, divShapeFunc] = simplexsf.uniformPShapeFunctionTable();
        %intergral error over the element 
        
        
    else
        string = [num2str(i),' is a transient element'];
        %disp(string)
        simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
        [ShapeFunc, divShapeFunc ] = simplexsf.variablePShapeFunctionTable();

    end
    
    u_numerical = zeros(1,length(IENall));
    for j = 1:length(IENall)
        u_numerical(j) = u(IENall(j));
    end
    
    uH = u_numerical*ShapeFunc;
    
    nP = size(qPoints,1);
    ux = zeros(nP,1);
    uy = zeros(nP,1);
    uquad = zeros(nP,1);
    
    tri = SimplexGeometry(IEN_mesh,vertexData);
    [JInverse, detJ,gij,xCord,yCord,A,hA] =  tri.jacobian();
    lam = [1-qPoints(:,1)-qPoints(:,2),qPoints(:,1:2)];
    
    for k = 1:nP
        uquad(k) = ShapeFunc(:,k)'*u_numerical';
        
        gradNGlobal = (divShapeFunc(:,:,k)*JInverse)';
        gradNx = gradNGlobal(1,:);
        gradNy = gradNGlobal(2,:);
        
        ux(k) = gradNx*u_numerical';
        uy(k) = gradNy*u_numerical';
    end
    
    
    x = xCord*lam';
    y = yCord*lam'; 
    [uA,ux_e,uy_e] = Exact(x,y,force,a,kappa);
    
    
    %Error = (((uH - uA)).^2*qPoints(:,3)*detJ)+Error;
    Error = (((uquad' - uA)).^2*qPoints(:,3)*detJ)+Error;
    H1_x = (((ux' - ux_e)).^2*qPoints(:,3)*detJ)+H1_x;
    H1_y = (((uy' - uy_e)).^2*qPoints(:,3)*detJ)+H1_y;
    %if Error>1e-6
        %string = [ num2str(i),' Order is ', num2str(pAll), ' error is ',num2str(Error^0.5) ];
        %disp(string)
    %end
        
end
Error = Error^0.5;
H1 = (H1_x+H1_y)^0.5;
end





