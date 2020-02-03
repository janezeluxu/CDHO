function [ ShapeFuncTable, divShapeFuncTable ] = ShapeFunc_quad(qPoints,p)
nQuadPoints = size(qPoints,1);
nShapeFunc = (p+1)*(p+1);
ShapeFuncTable = zeros(nShapeFunc,nQuadPoints);
divShapeFuncTable = zeros(nShapeFunc,2,nQuadPoints);
for i = 1: nQuadPoints
    IPoint = qPoints(i,:);
    [ ShapeFunc, divShapeFunc ] = QuaduniformP(IPoint,p,nShapeFunc);
    ShapeFuncTable(:,i) =  ShapeFunc;
    divShapeFuncTable(:,:,i) = divShapeFunc;
end
end

function [ ShapeFunc, divShapeFunc ] = QuaduniformP(IPoint,p,nShapeFunc)
xi = IPoint(1);
eta = IPoint(2);
        
[ SF1D_xi, diffSF1D_xi ] = ShapeFunc1D(xi,p);
[ SF1D_eta, diffSF1D_eta ] = ShapeFunc1D(eta,p);
ShapeFunc = zeros(nShapeFunc,1);
divShapeFunc = zeros(nShapeFunc,2);

% N1 to N16 for this quadrature point
%Order1Dto2D = [4,4;1,4;1,1;4,1;3,4;2,4;1,3;1,2;2,1;3,1;4,2;4,3;3,3;2,3;2,2;2,3];
%Order1Dto2D = [1,1;2,1;2,2;1,2];
Order1Dto2D = getOrder(p);
for i = 1:nShapeFunc
    order1 = Order1Dto2D(i,1);
    order2 = Order1Dto2D(i,2);
    ShapeFunc(i) = SF1D_xi(order1)*SF1D_eta(order2);
    divShapeFunc(i,1) = diffSF1D_xi(order1)*SF1D_eta(order2);
    divShapeFunc(i,2) = SF1D_xi(order1)*diffSF1D_eta(order2);
end
% ShapeFunc
% ShapeFunc = [(1-xi)*(1-eta);xi*(1-eta);xi*eta;(1-xi)*eta];
% divShapeFunc
% divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];
end

function Order1Dto2D = getOrder(p)
if p==1
Order1Dto2D = [1,1;2,1;2,2;1,2];
elseif p==2
    Order1Dto2D = [1,1;3,1;3,3;1,3;
                   2,1;3,2;2,3;1,2
                   2,2];
elseif p==3
    Order1Dto2D = [1,1;4,1;4,4;1,4;
        2,1;3,1;4,2;4,3;3,4;2,4;1,3;1,2;
        2,2;3,2;3,3;2,3];
elseif p==4
    Order1Dto2D = [1,1;5,1;5,5;1,5;
        2,1;3,1;4,1;5,2;5,3;5,4;4,5;3,5;2,5;1,4;1,3;1,2;
        2,2;3,2;4,2;2,3;3,3;4,3;2,4;3,4;4,4];
elseif p==5
    Order1Dto2D = [1,1;6,1;6,6;1,6;
        2,1;3,1;4,1;5,1
        6,2;6,3;6,4;6,5;
        5,6;4,6;3,6;2,6;
        1,5;1,4;1,3;1,2;
        2,2;3,2;4,2;5,2;
        2,3;3,3;4,3;5,3;
        2,4;3,4;4,4;5,4
        2,5;3,5;4,5;5,5];
elseif p==6
    Order1Dto2D = [1,1;7,1;7,7;1,7;
        2,1;3,1;4,1;5,1;6,1;
        7,2;7,3;7,4;7,5;7,6;
        6,7;5,7;4,7;3,7;2,7;
        1,6;1,5;1,4;1,3;1,2;
        2,2;3,2;4,2;5,2;6,2
        2,3;3,3;4,3;5,3;6,3
        2,4;3,4;4,4;5,4;6,4
        2,5;3,5;4,5;5,5;6,5
        2,6;3,6;4,6;5,6;6,6];
elseif p==7
    Order1Dto2D = [1,1;8,1;8,8;1,8;
        2,1;3,1;4,1;5,1;6,1;7,1;
        8,2;8,3;8,4;8,5;8,6;8,7;
        7,8;6,8;5,8;4,8;3,8;2,8;
        1,7;1,6;1,5;1,4;1,3;1,2;
        2,2;3,2;4,2;5,2;6,2;7,2;
        2,3;3,3;4,3;5,3;6,3;7,3
        2,4;3,4;4,4;5,4;6,4;7,4
        2,5;3,5;4,5;5,5;6,5;7,5;
        2,6;3,6;4,6;5,6;6,6;7,6
        2,7;3,7;4,7;5,7;6,7;7,7];
end

end

function [ShapeFunc,DivSF] = ShapeFunc1D(x,p)
L = length(x);
DivSF = zeros(p+1,L);
if p==1
    DivSF(1,:) = -1;
    DivSF(2,:) = 1;
else
    for i = 1:length(x)
        xi = x(i);
        Bpn1 = [0,bernsteinMatrix(p-1,xi),0];
        Bpn2 = [0,0,bernsteinMatrix(p-2,xi),0,0];
        dB = zeros(p+1,1);
        dB2 = zeros(p+1,1);
        for k = 1:p+1
            dB(k) = p*(Bpn1(k) - Bpn1(k+1));
            dB2(k) = p*(p-1)*(Bpn2(k)-2*Bpn2(k+1)+Bpn2(k+2));
        end
        DivSF(:,i) = dB;
    end
end
B = bernsteinMatrix(p,x);
ShapeFunc = B';
end