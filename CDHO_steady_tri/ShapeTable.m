function [ShapeFuncTable,divSFtable] = ShapeTable(nInt,p)

ShapeFuncTable = cell(length(p),1);  
divSFtable = cell(length(p),1);
for j = 1:length(p)
    n = nIntergerPoints(p(j),nInt);
    qPoints = TriGaussPoints(n);
    
    %[ ShapeFunc, divShapeFunc ] = ShapeFuncP2(qPoints);
    
    simplexsf = SimplexShapeFunc(qPoints,p(j));
    [ShapeFunc, divShapeFunc ] = simplexsf.uniformPShapeFunctionTable();  
    
    ShapeFuncTable{p(j)} =  ShapeFunc;
    divSFtable{p(j)} = divShapeFunc;
end
end