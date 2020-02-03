%clear all
%clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nsd; 
global kappa; 
%global direction
global force
%global BCType;
global a;
global Dirichlet_BC; 
global Dirichlet_edge;
global tol; 
global gamma;
global IntByPart;
%global StepFunc;
global itau; 
global stabflag;
global c2;
global gradMesh;
global nDirection;
global casenumber;
global static_factor;
global shapeFuncType;
global master;
global slave;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set problem variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = [1,0;0,1]*1e-4;
c2 = 1;
nDirection = [1;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsd = 2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solver settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-10; 
GaTol = 1e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IntByPart = 'no';
gamma = 1;
global CbI;
CbI = 4;

stabflag = 1;   % Add stabilization (1: VMS, 2: ...)
itau = 2; 
casenumber = 200;

if itau == 2
    itaucase ='static';
elseif itau == 0
    itaucase ='galerkin';
elseif itau == 3
    itaucase ='dynamic';
end
maxIter = 10;
%% mesh file
gradMesh = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if casenumber == 1
    force = 5;
    a = [0 0];
    c1_start = 0;
    Dirichlet_BC{1} = [2;18;14;16;5;6];
    Dirichlet_BC{2} = {0;0;0;0;0;0;};
    BCType = {'Strong','flux','Strong','flux'};
    meshFile =  "./mesh/s1.txt";
    casefile = string(casenumber);
    static_factor = 9;
    output_vtk = strcat('./testvtk/cases/',casefile,itaucase,'.vtk');
    Error = main(maxIter,meshFile,output_vtk,c1_start,GaTol);
elseif casenumber == 2
    %% case test variableP
    force = 6;
    a = [0 0];
    c1_start = 0;
    Dirichlet_BC{1} = [2;18;14;16;5;6];
    Dirichlet_BC{2} = {0;0;0;0;0;0;};
    BCType = {'Strong','flux','Strong','flux'};
    meshFile =  "./mesh/s2.txt";
    casefile = string(casenumber);
    static_factor = 9;
    output_vtk = strcat('./testvtk/cases/',casefile,itaucase,'.vtk');
    Error = main(maxIter,meshFile,output_vtk,c1_start,GaTol);
elseif casenumber == 3
    % case test variableP
    force = 7;
    a = [0 0];
    c1_start = 0;
    Dirichlet_BC{1} = [2;148;147;72;71;24;18;14;16;79;78;142;141;23;5;6];
    Dirichlet_BC{2} = {0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0};
    BCType = {'Strong','flux','Strong','flux'};
    meshFile =  "./mesh/s3.txt";
    casefile = string(casenumber);
    static_factor = 9;
    output_vtk = strcat('./testvtk/cases/',casefile,itaucase,'.vtk');
    Error = main(maxIter,meshFile,output_vtk,c1_start,GaTol);
elseif casenumber == 4
    % case StepFunc
    a = [1 1];
    c1_start = (1/norm(a));
    force = 0;
    Dirichlet_BC{1} = [21;93;101;2;6;15;5];
    Dirichlet_BC{2} = {0;0;0;1;1;0;1};
    order = 'p5';
    meshFile = strcat('./mesh/internal',order,'.txt')
    static_factor = 9;   
    casefile = "./internalLayer/";
    shapeFuncType = 'B';
    k = '4';
    output_filename = strcat('./testvtk/cases/',casefile,order,itaucase,'B.vtk');
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,k,'.txt')
    %main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    %plotSolution(meshFile,solutionFile)
    %getErrorwithRef(meshFile,solutionFile)
    figure(2)
    y = 0:0.01:1;
    x = ones(1,length(y))*0.2;
    LinePlot(meshFile,solutionFile,x,y)
    hold on
elseif casenumber == 5
    tic
    % boundary layer
    a =[1 0];
    %c1_start = (0.5/norm(a)/7);
    force = 0;
    Dirichlet_BC{1} = [184;190;15;16;10;5;6;2;117;110;185];
    Dirichlet_BC{2} = {0;0;1;1;1;1;1;1;0;0;0};
    %direction = 'X';
    %BCType = {'Strong','Strong','Strong','flux'};
    static_factor = 9;
    %meshFile = "./mesh/variableOrder666.txt";
    %% Dirichlet_BC for variable order mesh
    %Dirichlet_BC{1} = [184;190;15;16;10;5;6;2;117;110;331;330;400;399;185];
    %Dirichlet_BC{2} = {0;0;1;1;1;1;1;1;0;0;0;0;0;0;0};
    
    casefile = "./boundary/";
    order = 'p7';
    shapeFuncType = 'B';
    %order = 'adapt_p10';
    meshFile = strcat('./mesh/boundary',order,'.txt')
    k='4';
    output_filename = strcat('./testvtk/cases/',casefile,order,itaucase,'test.vtk');
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',casefile,order,itaucase,k,'.txt')
    main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    
    %getErrorwithRef(meshFile,solutionFile)
    %plotSolution(meshFile,solutionFile)
    
    %% plot at x=0.8
    figure(2)
    y = 0:0.002:0.2;
    x = ones(1,length(y))*0.82;
    LinePlot(meshFile,solutionFile,x,y)
    hold on
    toc
    
elseif casenumber == 6
    % case StepFunc
    a = [1 -1];
    c1_start = (1/norm(a));
    force = 0;
    Dirichlet_BC{1} = [21;93;101;2;6;5;11;10;16;15];
    Dirichlet_BC{2} = {0;0;0;1;1;1;1;1;0;0};
    order = 'p3ref1';
    meshFile = strcat('./mesh/internal',order,'.txt')
    static_factor = 9;   
    casefile = "./internalBL_corner/";
    shapeFuncType = 'B';
    k = '4';
    output_filename = strcat('./testvtk/cases/',casefile,order,itaucase,'B.vtk');
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,k,'.txt')
    main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    plotSolution(meshFile,solutionFile)
%     getErrorwithRef(meshFile,solutionFile)
%     figure(2)
%     x = 0:0.01:1;
%     y = ones(1,length(x))*0.2;
%     [concentration]=LinePlot(meshFile,solutionFile,x,y);
%     plot(x,concentration,'b*-')
%     hold on

elseif casenumber == 7
    % case StepFunc
    a = [1 1];
    c1_start = (1/norm(a));
    force = 0;
    Dirichlet_BC{1} = [79;83;35;36;30;25;26;20;134;130];
    Dirichlet_BC{2} = {0;0;1;1;1;0;0;0;0;0};
    order = 'p1';
    meshFile = strcat('./mesh/internalBL',order,'.txt')
    static_factor = 9;   
    casefile = "./internalBL/";
    shapeFuncType = 'B';
    k = '4';
    output_filename = strcat('./testvtk/cases/',casefile,order,itaucase,'B.vtk');
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,k,'.txt')
    main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    plotSolution(meshFile,solutionFile)
%     getErrorwithRef(meshFile,solutionFile)
%     figure(2)
%     x = 0:0.01:1;
%     y = ones(1,length(x))*0.2;
%     [concentration]=LinePlot(meshFile,solutionFile,x,y);
%     plot(x,concentration,'b*-')
%     hold on
elseif casenumber == 100
    % case verification
    
    a = [1 1]*(1/sqrt(2));
    force = 0;
    Dirichlet_BC{1} = [2;5;16;18;6;14];
    Dirichlet_BC{2} = {0;0;1;1;0;1};
    Dirichlet_edge = [6;14];
    static_factor = 9;
    master = [2;24;18];
    slave = [5;23;16];
    casefile = "./1Dinclined/";
    order = 'p72by2';
    meshFile = strcat('./mesh/caseverify',order,'.txt')
    output_filename = strcat('./testvtk/verifycase/',casefile,order,itaucase,'.vtk')
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,'.txt')
    main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    
    %getErrorwithExact(meshFile,solutionFile);
    plotnum = 2;
    plotSolution(meshFile,solutionFile,plotnum);
    
%     figure(2)
%     hold on
%     x = -0.2/sqrt(2):0.01:0.8/sqrt(2);
%     y = x+0.4/sqrt(2);
%     s = (x+y)/(sqrt(2));
%     concentration = LinePlot(meshFile,solutionFile,x,y);
%     Lineplotfilename = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,'line.txt')
%     writeSolution(concentration,Lineplotfilename)
%     plot(s,concentration,'r-o')

elseif casenumber == 200
    % case verification
    
    a = [1 0];
    force = 0;
    Dirichlet_BC{1} = [2;24;18;16;23;5];
    Dirichlet_BC{2} = {0;0;0;1;1;1};
    Dirichlet_edge = [24;23];
    static_factor = 9;
    master = [2;6;5];
    slave = [18;14;16];
    casefile = "./1D/";
    order = 'p72by2';
    meshFile = strcat('./mesh/caseverify1D',order,'.txt')
    output_filename = strcat('./testvtk/cases/',casefile,order,itaucase,'.vtk')
    solutionFile = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,'.txt')
    main(maxIter,meshFile,output_filename,solutionFile,GaTol)
    
    %getErrorwithExact(meshFile,solutionFile);
    plotnum = 3;
    plotSolution(meshFile,solutionFile,plotnum);
    
    %plotnum = 3;
    %exactSolution = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,'exact.txt')
    %getexactSolution(maxIter,meshFile,output_filename,exactSolution,GaTol);
    %plotSolution(meshFile,exactSolution,plotnum);
    figure(2)
    x = 0:0.001:1;
    y = ones(1,length(x))*0.2;
    concentration = LinePlot(meshFile,solutionFile,x,y);
    Lineplotfilename = strcat('./testvtk/cases/',casefile,'solution/',order,itaucase,'line.txt')
    writeSolution(concentration,Lineplotfilename)
    plot(x,concentration,'b-*')
    
end