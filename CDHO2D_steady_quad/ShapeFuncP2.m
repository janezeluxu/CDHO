function [ ShapeFunc, divShapeFunc ] = ShapeFuncP2(IntPoint)
lam2 = IntPoint(1);
lam3 = IntPoint(2);
lam1 = 1-lam2-lam3;

ShapeFunc = [lam1,lam2,lam3,-2*lam1*lam2,-2*lam2*lam3,-2*lam1*lam3];
divShapeFunc = [-1,-1;1,0;0,1;4*lam2 + 2*lam3 - 2,2*lam2; ...
    -2*lam3, -2*lam2;2*lam3,2*lam2 + 4*lam3 - 2;];

N1 = 2*(1-lam2-lam3)*(0.5-lam2-lam3);
N2 = 2*lam2*(lam2-0.5);
N3 = 2*lam3*(lam3-0.5);
N4 = 4*lam2*lam1;
N5 = 4*lam2*lam3;
N6 = 4*lam3*lam1;

dN1lam2 = 4*(lam2+lam3)-3;
dN2lam2 = 4*lam2 - 1;
dN3lam2 = 0;
dN4lam2 = 4 - 4*lam3 - 8*lam2;
dN5lam2 = 4*lam3;
dN6lam2 = -4*lam3;

dN1lam3 = 4*lam2 + 4*lam3 - 3;
dN2lam3 = 0;
dN3lam3 = 4*lam3 - 1;
dN4lam3 = -4*lam2;
dN5lam3 = 4*lam2;
dN6lam3 = 4 - 8*lam3 - 4*lam2;

ShapeFunc = [N1,N2,N3,N4,N5,N6];
divShapeFunc = [dN1lam2,dN1lam3;dN2lam2,dN2lam3;dN3lam2,dN3lam3;...
    dN4lam2,dN4lam3;dN5lam2,dN5lam3;dN6lam2,dN6lam3];
end