function polyplot(meshData,solution_ien,p_ien,edge_usage,vertexData,u,plotnum)
%global TotalDOF; 
global Grid_size;
%create uniform p shape function tables, p is a list
for ele = 1:Grid_size 
%for ele = 1:1     
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
    x4 = vertexData{vIDs(4),2}(1);
    y4 = vertexData{vIDs(4),2}(2);
    
    
    x = [x1,x2,x3,x4];
    y = [y1,y2,y3,y4];
            
    %% plot high order solution for each element 
    dx = 0:0.01:1;
    [xc,yc,uIntPoint,sx,sy,uexact] = ElementSurfacePlot(IENall,pAll,variable_ele,edgeU,dx,x,y);
    %subplot(2,2,plotnum)
    %sp = surf(xc,yc,uIntPoint);
    %set(sp,'EdgeColor',[1 1 1]);
    %set(sp,'FaceColor',[1 1 1]);
    %colormap default
    %sp.EdgeColor = 'none';
    
    sp = contour(xc,yc,uIntPoint,[0.001:0.001:0.05],'linewidth',2);
    %sp = contour(xc,yc,uexact,[0.1:0.2:1.0],'linewidth',2);
    %uIntPoint
    %sp = contour(sx,sy,uIntPoint,[0.1:0.2:1.0],'linewidth',2);
    %sp = contour(xc,yc,uIntPoint,[-0.2:0.04:1.2]);
    hold on
    
end
%light
%light('Position',[0.75 0.5 1.2],'Style','infinite')
end

function [xc,yc,uIntPoint,sx,sy,uexact] = ElementSurfacePlot(IENall,pAll,variable_ele,edgeU,dx,x,y)
global kappa; 
global casenumber;
global a;

[xi,eta] = meshgrid(dx,dx);
    
    uIntPoint = zeros(size(xi,1),size(eta,2));
    uexact = zeros(size(xi,1),size(eta,2));
    xc = zeros(size(xi,1),size(eta,2));
    yc = zeros(size(xi,1),size(eta,2));
    sx = zeros(size(xi,1),size(eta,2));
    sy = zeros(size(xi,1),size(eta,2));
    for i = 1:size(xi,1)
        for j = 1:size(eta,1)
            
            N1 = (1-xi(i,j))*(1-eta(i,j));
            N2 = xi(i,j)*(1-eta(i,j));
            N3 = xi(i,j)*eta(i,j);
            N4 = (1-xi(i,j))*eta(i,j);
            
            xc(i,j) = N1*x(1)+N2*x(2)+N3*x(3)+N4*x(4);
            yc(i,j) = N1*y(1)+N2*y(2)+N3*y(3)+N4*y(4);
    
            sx(i,j) = (xc(i,j)+yc(i,j))/(sqrt(2));
            sy(i,j) = abs(xc(i,j)-yc(i,j))/(sqrt(2));
            
            qPoints(1) = xi(i,j);
            qPoints(2) = eta(i,j);
            %simplexsf = SimplexShapeFunc(qPoints,p); %%assume uniform p
            %[ShapeFunc, ~ ] = simplexsf.uniformPShapeFunctionTable();  
            
            if range(pAll) == 0
                %disp('uniform p element');
                p=pAll(1);
                %simplexsf = SimplexShapeFunc(qPoints,p); %%assume uniform p
                %[ShapeFunc, ~ ] = simplexsf.uniformPShapeFunctionTable();  
                [ ShapeFunc, ~ ] = ShapeFunc_quad(qPoints,p);
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
            
            uexact(i,j) = Exact(xc(i,j),yc(i,j),casenumber,a,kappa);
        end
    end
end