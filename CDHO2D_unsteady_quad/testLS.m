
function testLS()
xLS =[0.707107000000000,   0.618718000000000,   0.618718000000000,   0.530330000000000,...
0.618718000000000   0.530330000000000   0.530330000000000   0.441942000000000]


yLS =[0.707107000000000   0.618718000000000   0.795495000000000   0.707107000000000,...
0.618718000000000   0.530330000000000   0.707107000000000   0.618718000000000]


uLS =[1.000000000000000   0.333584171218975   1.000000000000000   0.333587920973481,...
0.333584171218975   0.111211458891983   0.333587920973481   0.111212072163914]

s = (xLS+yLS)/(sqrt(2))
lhs = [ones(length(xLS),1),xLS',yLS'];
rhs = uLS';
lhs\rhs


lhs = [ones(length(xLS),1),s'];
rhs = uLS';
lhs\rhs
p = 1;
xmin = min(xLS); xmax=max(xLS); ymin=min(yLS);ymax=max(yLS);
[xi,eta]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax);
LS1 = getLS(xLS,yLS,uLS,p)
LS2 = getLS(xi,eta,uLS,p)

yy = [0,0,0.125,0.125,0,0,0.125,0.125];
xmin = min(s); xmax=max(s); ymin=min(yy);ymax=max(yy);
LS3 = getLS(s,yy,uLS,p)
[xi,eta]=mapxy(s,yy,xmin,xmax,ymin,ymax);
LS4 = getLS(xi,eta,uLS,p)

xeval =0.698814332603228;
yeval =0.707106955290478;


valuex_1 = LS1(1)+LS1(2)*xeval+LS1(3)*yeval
xmin = min(xLS); xmax=max(xLS); ymin=min(yLS);ymax=max(yLS);
[xi,eta]=mapxy(xeval,yeval,xmin,xmax,ymin,ymax);
valuex_2 = LS2(1)+LS2(2)*xi+LS2(3)*eta


sx = (xeval+yeval)/(sqrt(2));
sy = 0;
valuex_3 = LS3(1)+LS3(2)*sx+LS3(3)*sy

xmin = min(sx); xmax=max(sx); ymin=min(sy);ymax=max(sy);
[xi,eta]=mapxy(sx,sy,xmin,xmax,ymin,ymax);
valuex_4 = LS4(1)+LS4(2)*xi+LS4(3)*eta
end

function [xi,eta]=mapxy(xLS,yLS,xmin,xmax,ymin,ymax)
xi = (xLS-xmin)./(xmax-xmin);
eta = (yLS-ymin)./(ymax-ymin);
end

function LS = getLS(x,y,u,p)
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
end

