function [ue,ue_x,ue_y] = Exact(s,y,force,a,kappa)
global casenumber 
%x = 0:0.25:1;
%y = 0:0.25:1;
%x = (x+y)/(2^0.5)
%force=0;

ue = zeros(1,length(s));
if casenumber == 100
    g2 = 1;
    f = 0;
    a = norm(a);
    kappa = norm(kappa);
    r = -a/kappa;
    x = (s+y)/(sqrt(2));
    for i = 1:length(x)
        if x(i)>1          
            x(i)=1;
        end
    end
    [ue] = -(exp(r/2*(1-x))+exp(r/2))./(exp(r/2)+1).*(exp(r/2*(1-x))-exp(r/2))./(exp(r/2)-1);
    ue = (ue*g2+(f/a).*(x-ue));
    
    uex = (g2*r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) - exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) -...
        (f*((r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) + exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) +...
        (r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) - exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) - 1))/a ...
        + (g2*r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) + ...
        exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1));
    ue_x = uex/(sqrt(2));
    ue_y = uex/(sqrt(2));
    
elseif casenumber == 200
    g2 = 1;
    f = 0;
    a = norm(a);
    kappa = norm(kappa);
    r = -a/kappa;
    x = s;
    for i = 1:length(x)
        if x(i)>1          
            x(i)=1;
        end
    end
    [ue] = -(exp(r/2*(1-x))+exp(r/2))./(exp(r/2)+1).*(exp(r/2*(1-x))-exp(r/2))./(exp(r/2)-1);
    ue = (ue*g2+(f/a).*(x-ue));
    
    uex = (g2*r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) - exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) -...
        (f*((r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) + exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) +...
        (r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) - exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1)) - 1))/a ...
        + (g2*r*exp(-(r*(x - 1))./2).*(exp(-(r*(x - 1))/2) + ...
        exp(r/2)))/(2*(exp(r/2) - 1)*(exp(r/2) + 1));
    ue_x = uex;
    ue_y = 0;
elseif force ==1
    [ue] = (x.^2-1).*(y.^2-1);
elseif force ==2
    [ue] = (x.^2-1).*(y.^2-1);
elseif force == 3
    [ue] = SteadyHeatAnalytical(x,y);
elseif force ==4
    %uex = (1-x.^4);
    [ue] = sin(x).*sin(y);
elseif force ==5
    for i = 1:length(x)
        ue(i) = PieceWise1(x(i),y(i));
    end
elseif force ==6
    for i = 1:length(x)
        ue(i) = PieceWise2(x(i),y(i));
    end
elseif force ==7
    for i = 1:length(x)
        ue(i) = PieceWise3(x(i),y(i));
    end
elseif force ==8
    for i = 1:length(x)
        ue(i) = PieceWise4(x(i),y(i));
    end
end
end

function u = PieceWise1(x,y)
if x<= -0.6 
    u = (x+1);%*(1-y^2);%*(y+1);
elseif x>0.6 
    u = (1-x);%*(1-y^2);%*(1-y);
else
    u = (-(5/6)*x^2+0.7);%*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise2(x,y)
if x<= -0.6 
    u = (x+1);%*(1-y^2);%*(y+1);
elseif x>0.6 
    u = (1-x);%*(1-y^2);%*(1-y);
else
    a = -1/(4*0.6^3);
    b = 0.4-a*0.6^4;
    u = (a*x^4+b);%*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise3(x,y)
if x<= -0.6 
    u = (x+1)*(1-y^2);
elseif x>0.6 
    u = (1-x)*(1-y^2);
else
    u = (-(5/6)*x^2+0.7)*(1-y^2);
end
%u = (x^2-1);
end

function u = PieceWise4(x,y)
if (x> -0.6 && x< 0.6 && y> -0.6 && y< 0.6)
    u = (-x^2+0.6^2)*(0.6^2-y^2);
else
    u = 0;
end
end


