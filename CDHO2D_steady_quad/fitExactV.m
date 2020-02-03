function [vfitE] = fitExactV(vertexData,vID,IBC)
global a
SurroundNode = vertexData{vID,4};
vCordNode = vertexData{vID,2};
xc = vCordNode(1);
yc = vCordNode(2);
value = zeros(length(SurroundNode),1)-100;
for j = 1:length(SurroundNode)
    Snode = SurroundNode(j);
    if ismember(Snode,IBC)
        vCoodS = vertexData{Snode,2};
        xs = vCoodS(1);
        ys = vCoodS(2);
        D = [xs-xc,ys-yc];
        value(j) = D*a'/norm(D)/norm(a);
    end
end
[M,I] = max(value);
vfitE = SurroundNode(I);
end