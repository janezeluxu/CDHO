function plotHshape()
x = [0;1;0];
y = [0;0;1];
dx = 0:0.02:1;
p = 5;
for i = 1:(p+1)*(p+2)/2
    [xc,yc,uIntPoint] = ElementSurfacePlot(p,dx,x,y,i);
    sp = surf(xc,yc,uIntPoint);
    hold on
end
end

function [xc,yc,Na] = ElementSurfacePlot(p,dx,x,y,a)
    [t1,t2] = meshgrid(dx,dx);
    %[t1,t2] = meshgrid(qPoints(:,1),qPoints(:,2));
    lam1 = t1;
    lam2 = t2.*(1-t1);
    lam3 = (1-t1).*(1-t2);
    
    %N2 = lam2-0.5*2*lam2.*lam3;
    
    xc = lam1*x(1)+lam2*x(2)+lam3*x(3);
    yc = lam1*y(1)+lam2*y(2)+lam3*y(3);
    Na = zeros(size(lam2,1),size(lam2,2));
    for i = 1:size(lam2,1)
        for j = 1:size(lam2,2)
            qPoints(1) = lam2(i,j);
            qPoints(2) = lam3(i,j);
            [ ShapeFunc, divShapeFunc ] = HierarchicalShapeFunc(p,qPoints);
            Na(i,j) = ShapeFunc(a);
        end
    end
end