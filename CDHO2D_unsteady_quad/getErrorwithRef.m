
function getErrorwithRef(fileNameM1,solutionFileM1)
casefile = "./internalBC/";
order = 'p10ref3';
k = '4';
itaucase = 'dynamic';
fileNameM2 = strcat('./mesh/internal',order,'.txt')
solutionFileM2 = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,k,'.txt')
    
% casefile = "./internalLayer/";
% k = '3';
% order = 'p10ref3';
% itaucase = 'dynamic';
% fileNameM2 = strcat('./mesh/internal',order,'.txt')
% solutionFileM2 = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,k,'.txt')
%     
% casefile = "./boundary/";
% order = 'p10ref3';
% fileNameM2 = strcat('./mesh/boundary',order,'.txt')
% k='4';
% solutionFileM2 = strcat('./testvtk/cases/',casefile,'solution/',casefile,order,'dynamic',k,'.txt');

[L2Error,H1Error]=...
    getError(fileNameM1,fileNameM2,solutionFileM1,solutionFileM2)

end


function [L2Error,H1Error]=getError(fileNameM1,fileNameM2,solutionFileM1,solutionFileM2)
refineLevel = 6;
nInt = 20;
global Grid_size;
    
%% get error with reference solution. map the coarse mesh solution to fine 
% mesh. M1 is the coarse mesh, M2 is the fine mesh, every M1 mesh consist
% of 2^refineLevel numbers of M2 mesh element

[pM2,meshDataM2,vertexDataM2,Area,IBC,IBC_v,BCval,BCb,uHBCE,solution_ienM2,...
    p_ienM2,iper] = preprocess(fileNameM2);
Grid_size
[pM1,meshDataM1,vertexDataM1,Area,IBC,IBC_v,BCval,BCb,uHBCE,solution_ienM1,...
    p_ienM1,iper] = preprocess(fileNameM1);

uM1 = load(solutionFileM1);
uM2 = load(solutionFileM2);
Grid_sizeM1 = Grid_size
[L2Error,H1Errorx,H1Errory]=calcuError(pM1,pM2,uM1,uM2,nInt,Grid_sizeM1,refineLevel,...
                meshDataM2,solution_ienM2,p_ienM2,vertexDataM2,...
                meshDataM1,solution_ienM1,p_ienM1,vertexDataM1);
H1Error =  ( H1Errorx+H1Errory )^0.5;  
L2Error = L2Error^0.5;
end

function [L2Error,H1Errorx,H1Errory]=calcuError(pM1,pM2,uM1,uM2,nInt,Grid_sizeM1,refineLevel,...
                meshDataM2,solution_ienM2,p_ienM2,vertexDataM2,...
                meshDataM1,solution_ienM1,p_ienM1,vertexDataM1)

%% get error with reference solution. map the coarse mesh solution to fine mesh

Nelement = 2^refineLevel;
[ShapeFuncTableM2,divSFtableM2] = ShapeTable(nInt,pM2);
L2Error = 0;
H1Errorx = 0;
H1Errory = 0;
for eleM1 = 1:Grid_sizeM1
    %eleM1
    if mod(eleM1,32)==0
        eleM1
    end
    [IEN_meshM1,IENallM1,pAllM1] = elementIEN(eleM1,meshDataM1,solution_ienM1,p_ienM1);
    tri = SimplexGeometry(IEN_meshM1,vertexDataM1);
    [JInverseM1, detJM1,~,xCordM1,yCordM1,~,~] =  tri.jacobian();
    
    uEM1 = zeros(1,length(IENallM1));
    for j = 1:length(IENallM1)
        uEM1(j) = uM1(IENallM1(j));
    end
    H1(eleM1) = 0;  
    for ni = 1:Nelement
        eleM2 = (eleM1-1)*Nelement+ni;
        [IEN_mesh,IENall,pAll] = elementIEN(eleM2,meshDataM2,solution_ienM2,p_ienM2);
        nssl = length(IENall);
        n = nIntergerPoints(max(pAll),nInt);
        qPoints = TriGaussPoints(n);
        
        if range(pAll) == 0
            %disp('uniform p element');
            ShapeFunc = ShapeFuncTableM2{pAll(1)};
            divShapeFunc = divSFtableM2{pAll(1)};
            
            simplexsf = SimplexShapeFunc(qPoints,pAll(1));
            [ShapeFunc, divShapeFunc] = simplexsf.uniformPShapeFunctionTable(); 
        else
            %disp('nonuniform p element');
            simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
            [ShapeFunc, divShapeFunc] = simplexsf.variablePShapeFunctionTable();            
        end
        
        uEM2 = zeros(1,length(IENall));
        for j = 1:length(IENall)
            uEM2(j) = uM2(IENall(j));
        end
        tri = SimplexGeometry(IEN_mesh,vertexDataM2);
        [JInverse, detJ,gij,xCord,yCord,A,hA] =  tri.jacobian();
            
        % get the list of u at quadrature points
        nP = size(qPoints,1); 
        uquadM2 = zeros(nP,1);
        duquadM2 = zeros(nP,2);
        for k = 1:nP
            gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
            NaGlobal = ShapeFunc(:,k);
            uquadM2(k) = uEM2*NaGlobal;
            duquadM2(k,:) = uEM2*gradNaGlobal;
        end
        
        %% get uM1 and duM1 at the quadrature points of M2,
        % first map quadrature points in M2 to x y coordinates
        % then to points in area quadinate in M1
        lam = [1-qPoints(:,1)-qPoints(:,2),qPoints(:,1:2)];
        x = xCord*lam';
        y = yCord*lam';
        evaluatePoints = getPoints(x,y,xCordM1,yCordM1); % map to area quadinates
        % get shape function values at evaluatePoints        
        if range(pAllM1) == 0
            %disp('uniform p element');
            simplexsf = SimplexShapeFunc(evaluatePoints,pAllM1(1));
            [ShapeFuncM1, divShapeFuncM1] = simplexsf.uniformPShapeFunctionTable();  
        else
            %disp('nonuniform p element');
            simplexsf = SimplexShapeFunc(evaluatePoints,IENallM1,pAllM1);
            [ShapeFuncM1, divShapeFuncM1] = simplexsf.variablePShapeFunctionTable();            
        end
        
        nP = size(qPoints,1);  
        uquadM1 = zeros(nP,1);
        duquadM1 = zeros(nP,2);
        
        for k = 1:nP
            gradNaGlobal = divShapeFuncM1(:,:,k)*JInverseM1;
            NaGlobal = ShapeFuncM1(:,k);
            uquadM1(k) = uEM1*NaGlobal;
            duquadM1(k,:) = uEM1*gradNaGlobal;
        end
        
        L2Error = ((uquadM1' - uquadM2').^2*qPoints(:,3)*detJ)+L2Error;
        H1Errorx = ((duquadM1(:,1)' - duquadM2(:,1)').^2*qPoints(:,3)*detJ)+H1Errorx;
        H1Errory = ((duquadM1(:,2)' - duquadM2(:,2)').^2*qPoints(:,3)*detJ)+H1Errory;
    
        H1(eleM1) = ((duquadM1(:,1)' - duquadM2(:,1)').^2*qPoints(:,3)*detJ)+...
            ((duquadM1(:,2)' - duquadM2(:,2)').^2*qPoints(:,3)*detJ)+H1(eleM1);
    end
end
end

function evaluatePoints=getPoints(x,y,xCordM1,yCordM1)
%map to area quadinates
evaluatePoints = zeros(length(x),2);
for i = 1:length(x)
    x1 = xCordM1(1);
    x2 = xCordM1(2);
    x3 = xCordM1(3);
    y1 = yCordM1(1);
    y2 = yCordM1(2);
    y3 = yCordM1(3);
    Atotal = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = x(i);
    x2 = xCordM1(2);
    x3 = xCordM1(3);
    y1 = y(i);
    y2 = yCordM1(2);
    y3 = yCordM1(3);
    A1 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = xCordM1(1);
    x2 = x(i);
    x3 = xCordM1(3);
    y1 = yCordM1(1);
    y2 = y(i);
    y3 = yCordM1(3);
    A2 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    x1 = xCordM1(1);
    x2 = xCordM1(2);
    x3 = x(i);
    y1 = yCordM1(1);
    y2 = yCordM1(2);
    y3 = y(i);
    A3 = (0.5*(-(x3-x1)*(y2-y1)+(x2-x1)*(y3-y1)));
    lam2 = A2/Atotal;
    lam3 = A3/Atotal;
    evaluatePoints(i,:) = [lam2,lam3];
end
%% check if A1/Atotal = 1-lam2-lam3
A1test = A1/Atotal;
lam1 = 1-lam2-lam3;
end