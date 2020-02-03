function [ShapeFuncTable,divSFtable] = ShapeTable(nInt,p)

    n = nIntergerPoints(p,nInt);
    qPoints = QuadGaussPoints(n);
    
    [ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(qPoints,p);
    
end