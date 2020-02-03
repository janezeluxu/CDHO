global stabflag; 
            global kappa;
            global a; 
            global force; 
            global IntByPart;
            global itau;
            global nsd; 
 
stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 2;
a = [1 -1];
c1_start = (1/norm(a));
force = 0;

kappa = [1,0;0,1]*1e-4;
c2 = 1;
nDirection = [1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-10;
GaTol = 5e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IntByPart = 'no';

ele = 33;
        IENall = [20    21    26];
        IEN_mesh = [20    21    26];
        pAll = [1 1 1 1];
        n = 2;
        qPoints = 0;
        vertexData = 0;
        
        p=1;
        [ShapeFuncTable,divSFtable] = ShapeTable(0,p);
        ShapeFunc = ShapeFuncTable{1};
        divShapeFunc = divSFtable{1};
        tauELE = 1;
        
        
        simplex = SimplexStiffMatrix(ele,IEN_mesh,IENall,pAll,n,qPoints,vertexData,...
        ShapeFunc,divShapeFunc,tauELE);
    %ShapeFunc
    %divShapeFunc
    [elementK,elementF] = simplex.eleStiffMatrix()