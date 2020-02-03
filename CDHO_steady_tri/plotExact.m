function plotExact()
[a,b] = meshgrid(0:0.1:1,0:0.1:1);
x = (a-b)/sqrt(2);
y = (a+b)/sqrt(2);
force = 0;
a = [1 1]*(1/sqrt(2));
kappa = [1,0;0,1]*1e-3;

%for i = 1:length()
[ue,ue_x,ue_y] = Exact(x,y,force,a,kappa)
surf(x,y,ue)
end