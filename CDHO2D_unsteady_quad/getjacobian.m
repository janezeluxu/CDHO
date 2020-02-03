function [JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
    getjacobian(ien_mesh,vertexData,xi,eta)
vIDs = ien_mesh;

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

hA = sqrt((x1-x2)^2+(y1-y2)^2);
A = hA^2;

J = [0,0;0,0];
dNl1 = [-(1-eta), 1-eta, eta,-eta]; % derivative of N1, N2 and N4 with xi
dNl2 = [-(1-xi), -xi, xi,(1-xi)]; % derivative of N1, N2 and N4 with eta

%divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];
for i = 1:4
    J(1,1) = J(1,1)+xCord(i)*dNl1(i); % derivative of x with xi
    J(1,2) = J(1,2)+xCord(i)*dNl2(i); % derivative of x with eta
    J(2,1) = J(2,1)+yCord(i)*dNl1(i); % derivative of y with xi
    J(2,2) = J(2,2)+ yCord(i)*dNl2(i); % derivative of y with eta
end
JInverse = J^-1;
gij = JInverse'*JInverse;
% Jp11 = JInverse(1,1);
% Jp12 = JInverse(1,2);
% Jp21 = JInverse(2,1);
% Jp22 = JInverse(2,2);
% 
% gij(1,1) =(Jp11)*(Jp11)+(Jp21)*(Jp21);
% gij(1,2) =(Jp11)*(Jp12)+(Jp21)*(Jp22);
% gij(2,1) =(Jp12)*(Jp11)+(Jp22)*(Jp21);
% gij(2,2) =(Jp12)*(Jp12)+(Jp22)*(Jp22);
detJ = det(J);
end