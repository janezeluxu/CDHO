function u = PieceWise1(x,y)
if x<= -0.6 
    u = (x+1);
elseif x>0.6 
    u = (1-x);
else
    u = (-(5/6)*x^2+0.7);
end
%u = (x^2-1);
end

function u = PieceWise2(x,y)
if x<= -0.6 
    u = (x+1);
elseif x>0.6 
    u = (1-x);
else
    a = -1/(4*0.6^3);
    b = 0.4-a*0.6^4;
    u = (a*x^4+b);
end
end

function u = PieceWise3(x,y)
if x<= -0.6 
    u = (x+1)*(1-y^2);
elseif x>0.6 
    u = (1-x)*(1-y^2);
else
    u = (-(5/6)*x^2+0.7)*(1-y^2);
end
end

function u = PieceWise4(x,y)
if (x> -0.6 && x< 0.6 && y> -0.6 && y< 0.6)
    u = (-x^2+0.6^2)*(0.6^2-y^2);
else
    u = 0;
end
end