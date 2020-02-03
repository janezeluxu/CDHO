function twod_to_vtk_tau (output_filename )
  timestamp ( );
  fprintf ( '\n' );
  fprintf ( 1, 'TWOD_TO_VTK_TEST:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Test the TWOD_TO_VTK library.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Read in TWOD finite element data.\n' );
  fprintf ( 1, '  Write it out as a legacy VTK file.\n' );
%
%  Load the data from text files.
%
  nodes = load ( './testvtk/nodes.txt' );
  elements = load ( './testvtk/elements.txt' );
  taus = load ( './testvtk/taus.txt' );
  tauDC = load ( './testvtk/tauDC.txt' );
  taut_percent = load ( './testvtk/taut_percent.txt' );
  tauadv_percent = load ( './testvtk/tauadv_percent.txt' );
  taudiff_percent = load ( './testvtk/taudiff_percent.txt' );
 %output_filename = './testvtk/verifycase/3.vtk';

  title = 'UVP for 1D advection-diffusion problem';
%
%  Have the data written to a file.
%
  twod_to_vtk_cell ( nodes, elements, taus,tauDC,taut_percent,tauadv_percent, taudiff_percent,output_filename, title );
%
%  Terminate.
%
  fprintf ( '\n' );
  fprintf ( 1, 'TWOD_TO_VTK_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( '\n' );
  timestamp ( );

  return
end

function twod_to_vtk_cell ( xy, e_conn, tau1,tau2,tau3,tau4, tau5, output_filename, title )
  [ node_num, dim_num ] = size ( xy );
  [ element_num,  element_order ] = size ( e_conn );
%
%  Open the output file.
%
  if ( isempty ( output_filename ) )
    output_filename = 'ns2d_fem.vtk';
  end

  output_unit = fopen ( output_filename, 'w' );
%
%  Transpose or otherwise modify the data.
%
  xyz = zeros ( 3, node_num );
  xyz(1:2,1:node_num) = xy(1:node_num,1:2)';

  if ( element_order == 6 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'TWO_TO_VTK:\n' );
    fprintf ( 1, '  The input data uses quadratic elements.\n' );
    fprintf ( 1, '  The output data will use linear elements.\n' );
  end

  element_order2 = 4;

  element_node = zeros ( element_order2, element_num );
  element_node(1:element_order2,1:element_num) = ...
    e_conn(1:element_num,1:element_order2)' - 1;

  taus = zeros ( 1, element_num );
  taus(1,1:element_num) = tau1(1:element_num,1)';
  
  tauDC = zeros ( 1, element_num );
  tauDC(1,1:element_num) = tau2(1:element_num,1)';

  taut_percent = zeros ( 1, element_num );
  taut_percent(1,1:element_num) = tau3(1:element_num,1)';
  
  tauadv_percent = zeros ( 1, element_num );
  tauadv_percent(1,1:element_num) = tau4(1:element_num,1)';
  
  taudiff_percent = zeros ( 1, element_num );
  taudiff_percent(1,1:element_num) = tau5(1:element_num,1)';
%
%  Write the data.
%
  vtk_puv_write ( output_unit, title, node_num, element_num, ...
    element_order2, xyz, element_node, taus,tauDC,taut_percent,tauadv_percent, taudiff_percent);

  fclose ( output_unit );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  The data was written to "%s"\n', output_filename );

  return
end

function vtk_puv_write ( output_unit, title, node_num, element_num, ...
  element_order, xyz, element_node, taus,tauDC,taut_percent,tauadv_percent, taudiff_percent)
  fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
  fprintf ( output_unit, '%s\n', title );
  fprintf ( output_unit, 'ASCII\n' );
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
  fprintf ( output_unit, 'POINTS %d double\n', node_num );

  for node = 1 : node_num
    fprintf ( output_unit, '  %2.16f  %2.16f  %2.16f\n', xyz(1:3,node) );
  end
%
%  Note that CELL_SIZE uses ELEMENT_ORDER+1 because the order of each element
%  is included as a data item.
%
  cell_size = element_num * ( element_order + 1 );

  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
  for element = 1 : element_num
    fprintf ( output_unit, '  %d', element_order );
    for order = 1 : element_order
      fprintf ( output_unit, '  %d', element_node(order,element) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  VTK has a cell type 22 for quadratic triangles.  However, we
%  are going to strip the data down to linear triangles for now,
%  which is cell type 5.
%
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELL_TYPES %d\n', element_num );

  if ( element_order == 3 )
    for element = 1 : element_num
      fprintf ( output_unit, '5\n' );
    end
  elseif ( element_order == 6 )
    for element = 1 : element_num
      fprintf ( output_unit, '22\n' );
    end
  elseif ( element_order == 4 )
    for element = 1 : element_num
      fprintf ( output_unit, '9\n' );
    end
  end

%   fprintf ( output_unit, '\n' );
%   fprintf ( output_unit, 'POINT_DATA %d\n', node_num );
%   fprintf ( output_unit, 'SCALARS concentration double\n' );
%   fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
%   for node = 1 : node_num
%     fprintf ( output_unit, '  %2.16f\n', concentration(1,node) );
%   end
%   
%   fprintf ( output_unit, 'SCALARS ga double\n' );
%   fprintf ( output_unit, 'LOOKUP_TABLE default\n' ); 
%   for node = 1 : node_num
%     fprintf ( output_unit, '  %2.16f\n', ga(1,node) );
%   end
%   fprintf ( output_unit, 'SCALARS tau double\n' );
%   fprintf ( output_unit, 'LOOKUP_TABLE default\n' ); 
%   for node = 1 : node_num
%     fprintf ( output_unit, '  %2.16f\n', tau(1,node) );
%   end
%   
% %   fprintf ( output_unit, 'VECTORS velocity double\n' );
% %   for node = 1 : node_num
% %     fprintf ( output_unit, '  %f  %f  %f\n', uvw(1:3,node) );
% %   end

fprintf ( output_unit, '\n' );
fprintf ( output_unit, 'CELL_DATA %d\n', element_num );
fprintf ( output_unit, 'SCALARS taus double\n' );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for node = 1 : element_num
    fprintf ( output_unit, '  %2.16f\n', taus(1,node) );
end

fprintf ( output_unit, 'SCALARS tauDC double\n' );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for node = 1 : element_num
    fprintf ( output_unit, '  %2.16f\n', tauDC(1,node) );
end

fprintf ( output_unit, 'SCALARS taut_percent double\n' );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for node = 1 : element_num
    fprintf ( output_unit, '  %2.16f\n', taut_percent(1,node) );
end

fprintf ( output_unit, 'SCALARS tauadv_percent double\n' );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for node = 1 : element_num
    fprintf ( output_unit, '  %2.16f\n', tauadv_percent(1,node) );
end

fprintf ( output_unit, 'SCALARS taudiff_percent double\n' );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for node = 1 : element_num
    fprintf ( output_unit, '  %2.16f\n', taudiff_percent(1,node) );
end
  
  return
end
