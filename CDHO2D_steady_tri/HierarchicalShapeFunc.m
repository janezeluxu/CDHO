function [ ShapeFunc, divShapeFunc ] = HierarchicalShapeFunc(p,IntPoint)
[ ShapeFunc, divShapeFunc ] = Hierarchical(p,IntPoint);
if p==2
   [ ShapeFunc, divShapeFunc ] = HierarchicalP2(IntPoint);
end
end
function [ ShapeFunc, divShapeFunc ] = HierarchicalP2(IntPoint)
lam2 = IntPoint(1);
lam3 = IntPoint(2);
lam1 = 1-lam2-lam3;

ShapeFunc = [lam1,lam2,lam3,-2*lam1*lam2,-2*lam2*lam3,-2*lam1*lam3];
divShapeFunc = [-1,-1;1,0;0,1;4*lam2 + 2*lam3 - 2,2*lam2; ...
    -2*lam3, -2*lam2;2*lam3,2*lam2 + 4*lam3 - 2;];
end

function [ ShapeFunc, divShapeFunc ] = Hierarchical(p,IntPoint)
%give one shape function and derivative value at one intergral point
%and one lexico order combination in lambda space
lam2 = IntPoint(1);
lam3 = IntPoint(2);
lam1 = 1-lam2-lam3;
lam = [lam1,lam2,lam3];
dlam = [-1,-1;1,0;0,1];

nsSF = (p+1)*(p+2)/2;
ShapeFunc = zeros(nsSF,1);
divShapeFunc = zeros(nsSF,2);

%% vertex mode
for i = 1: 3
    ShapeFunc(i) = lam(i);
    divShapeFunc(i,1) = dlam(i,1);
    divShapeFunc(i,2) = dlam(i,2);
end

%% edge mode 
count = 4;
k=0;
for order = 2:p
    [SF,divSF]=EdgeHierarchical(order,lam,dlam);
    for i = 1:3
        index = 4+(p-1)*(i-1)+k;
        ShapeFunc(index) = SF(i);
        divShapeFunc(index,1) = divSF(i,1);
        divShapeFunc(index,2) = divSF(i,2);
        count= count+1;
    end
    k=k+1;
end
%% face mode
for order=3:p
    [SF,divSF] = FaceHierarchical(order,lam);
    num = length(SF);
    for i = 1:num
        ShapeFunc(count) = SF(i);
        divShapeFunc(count,1) = divSF(i,1);
        divShapeFunc(count,2) = divSF(i,2);
        count= count+1;
    end
end
end

function [SF,divSF]=EdgeHierarchical(p,lam,dlam)
[fai1,zai1,dfai1,dzai1]=EdgeFunction(p,lam(1),lam(2),dlam(1,:),dlam(2,:));
[fai2,zai2,dfai2,dzai2]=EdgeFunction(p,lam(2),lam(3),dlam(2,:),dlam(3,:));
[fai3,zai3,dfai3,dzai3]=EdgeFunction(p,lam(3),lam(1),dlam(3,:),dlam(1,:));
SF = [fai1*zai1;fai2*zai2;fai3*zai3];
divSF = [dfai1*zai1+fai1*dzai1;dfai2*zai2+fai2*dzai2;dfai3*zai3+fai3*dzai3];
end

function [fai,zai,dfai,dzai]=EdgeFunction(p,lam1,lam2,dlam1,dlam2)
if p==2
    fai = 1;
    dfai = [0,0];
elseif p==3
    fai = (lam2-lam1);
    dfai = dlam2-dlam1;
elseif p==4
    fai = (lam1^2-3*lam1*lam2+lam2^2);
    dfai = 2*lam1*dlam1-3*lam1*dlam2-3*dlam1*lam2+2*lam2*dlam2;
elseif p==5
    fai = (lam2^3-6*lam1*lam2^2+6*lam1^2*lam2-lam1^3);
    dfai = 3*lam2^2*dlam2-6*dlam1*lam2^2-6*lam1*2*lam2*dlam2+6*2*lam1*dlam1*lam2+6*lam1^2*dlam2-3*lam1^2*dlam1;
end
zai = -2*lam1*lam2;
dzai = -2*lam1*dlam2-2*lam2*dlam1;
end

function [SF,divSF]=FaceHierarchical(p,lam)
lam1 = lam(1);
lam2 = lam(2);
lam3 = lam(3);
if p==3
    SF = lam1*lam2*lam3;
    divSF = [- lam3*(lam2 + lam3 - 1) - lam2*lam3,...
        - lam2*(lam2 + lam3 - 1) - lam2*lam3];
elseif p==4
    SF = lam1*lam2*lam3*[lam1-1/3;lam2-1/3];
    divSF = [(lam3*(9*lam2^2 + 12*lam2*lam3 - 10*lam2 + 3*lam3^2 - 5*lam3 + 2))/3,...
             (lam2*(3*lam2^2 + 12*lam2*lam3 - 5*lam2 + 9*lam3^2 - 10*lam3 + 2))/3;...
               -(lam3*(6*lam2*lam3 - lam3 - 8*lam2 + 9*lam2^2 + 1))/3,...
               -(lam2*(3*lam2 - 1)*(lam2 + 2*lam3 - 1))/3;];
elseif p==5
    SF = lam1*lam2*lam3*[lam1^2-(3/4)*lam1+(3/28);...
        lam1*lam2-(1/4)*(lam1+lam2)+(1/14); lam2^2-(3/4)*lam2+(3/28)];
    divSF = [-(lam3*(112*lam2^3 + 252*lam2^2*lam3 - 189*lam2^2 + 168*lam2*lam3^2 - 252*lam2*lam3 + 90*lam2 + 28*lam3^3 - 63*lam3^2 + 45*lam3 - 10))/28,...
              -(lam2*(28*lam2^3 + 168*lam2^2*lam3 - 63*lam2^2 + 252*lam2*lam3^2 - 252*lam2*lam3 + 45*lam2 + 112*lam3^3 - 189*lam3^2 + 90*lam3 - 10))/28;...
              (lam3*(112*lam2^3 + 168*lam2^2*lam3 - 168*lam2^2 + 56*lam2*lam3^2 - 126*lam2*lam3 + 66*lam2 - 7*lam3^2 + 12*lam3 - 5))/28,...
              (lam2*(28*lam2^3 + 112*lam2^2*lam3 - 56*lam2^2 + 84*lam2*lam3^2 - 126*lam2*lam3 + 33*lam2 - 21*lam3^2 + 24*lam3 - 5))/28;...
              -(lam3*(48*lam2 + 3*lam3 - 42*lam2*lam3 + 84*lam2^2*lam3 - 147*lam2^2 + 112*lam2^3 - 3))/28,...
              -(lam2*(28*lam2^2 - 21*lam2 + 3)*(lam2 + 2*lam3 - 1))/28;];
end
end