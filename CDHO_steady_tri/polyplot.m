function polyplot(meshData,solution_ien,p_ien,edge_usage,vertexData,u)
%global TotalDOF; 
global Grid_size;
%create uniform p shape function tables, p is a list
for ele = 1:Grid_size 
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions  
    %% get meshIEN, solutionIEN, and shape functions
    [IEN_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    nssl = length(IENall);
    edgeU = edge_usage(ele,:);
    variable_ele = zeros(nssl,1);
    
    for i = 1:nssl
        variable_ele(i) = u(IENall(i));
    end
    vIDs = IEN_mesh;
    x1 = vertexData{vIDs(1),2}(1);
    y1 = vertexData{vIDs(1),2}(2);
    x2 = vertexData{vIDs(2),2}(1);
    y2 = vertexData{vIDs(2),2}(2);
    x3 = vertexData{vIDs(3),2}(1);
    y3 = vertexData{vIDs(3),2}(2);
    
    
    x = [x1,x2,x3];
    y = [y1,y2,y3];
            
    %% plot high order solution for each element 
    dx = 0:0.2:1;
    [xc,yc,uIntPoint] = ElementSurfacePlot(IENall,pAll,variable_ele,edgeU,dx,x,y);
    %sp = surf(xc,yc,uIntPoint);
    %set(sp,'FaceColor',[1 1 1]);
    %sp.EdgeColor = 'none';
    
    %sp = contour(xc,yc,uIntPoint,[0.1:0.1:1.2],'linewidth',2);
    sp = contour(xc,yc,uIntPoint,[-0.2:0.08:0.9,1.1:0.08:1.2],'linewidth',2);
    hold on
    
end
%light
%light('Position',[0.75 0.5 1.2],'Style','infinite')
end

function [xc,yc,uIntPoint] = ElementSurfacePlot(IENall,pAll,variable_ele,edgeU,dx,x,y)
    [t1,t2] = meshgrid(dx,dx);
    %[t1,t2] = meshgrid(qPoints(:,1),qPoints(:,2));
    lam1 = t1;
    lam2 = t2.*(1-t1);
    lam3 = (1-t1).*(1-t2);
    
    %N2 = lam2-0.5*2*lam2.*lam3;
    
    xc = lam1*x(1)+lam2*x(2)+lam3*x(3);
    yc = lam1*y(1)+lam2*y(2)+lam3*y(3);
    uIntPoint = zeros(size(lam2,1),size(lam2,2));
    for i = 1:size(lam2,1)
        for j = 1:size(lam2,2)
            qPoints(1) = lam2(i,j);
            qPoints(2) = lam3(i,j);
            %simplexsf = SimplexShapeFunc(qPoints,p); %%assume uniform p
            %[ShapeFunc, ~ ] = simplexsf.uniformPShapeFunctionTable();  
            
            if range(pAll) == 0
                %disp('uniform p element');
                p=pAll(1);
                simplexsf = SimplexShapeFunc(qPoints,p); %%assume uniform p
                [ShapeFunc, ~ ] = simplexsf.uniformPShapeFunctionTable();  
            else
                %disp('nonuniform p element');
                simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
                [ShapeFunc, ~] = simplexsf.variablePShapeFunctionTable();
                
            end
    
            count = 3;
            for k = 1:3
                for l = 1:p-1
                    count = count+1;
                    ShapeFunc(count,:)= edgeU(k)*ShapeFunc(count,:);
                end
            end
            uIntPoint(i,j) = variable_ele'*ShapeFunc;
        end
    end
end