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
        if isa(a{1}, 'function_handle')==1
            velocity1 = a{1}(xs,ys);
        else
            velocity1 = a{1};
        end
        if isa(a{2}, 'function_handle')==1
            velocity2 = a{2}(xs,ys);
        else
            velocity2 = a{2};
        end
        velocity = [velocity1,velocity2];
        
        if norm(velocity)>0
        value(j) = D*velocity'/norm(D)/norm(velocity);
        else
            value(j) = 0;
        end
    end
end
[M,I] = max(value);
vfitE = SurroundNode(I);
end