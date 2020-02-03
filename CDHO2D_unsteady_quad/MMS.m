function [Source] = MMS(t,x,y,casenumber)
if casenumber == 300
    x = 10*x-1;
    y = 10*y-1;
    if (sqrt(x^2+y^2)>1)
        Source = 0;
    else
        Source = exp(-t^1)*cos((pi/2)*sqrt(x^2+y^2));
    end
elseif casenumber == 400
    Source = 0;
end

end