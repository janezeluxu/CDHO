function [qPoints,xi] = QuadGaussPoints(n)
[xi, w] = GaussQuad(n,1);
lxi = length(xi);
qPoints = zeros(lxi*lxi,3);
count = 1;
for i = 1:lxi
    for m = 1:lxi
        %n = (i-1)*lxi+m;
        qPoints(count,:)=[xi(i),xi(m),w(i)*w(m)];
        count = count+1;
    end
end

end