function twod_to_vtk_test (output_filename )

%*****************************************************************************80
%
%% TWOD_TO_VTK_TEST tests the TWOD_TO_VTK library.
%
%  Discussion:
%
%    The data represents a flow solution in a 10 by 3 channel.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 December 2010
%
%  Author:
%
%    John Burkardt
%
  global static_factor;
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
  uvp = load ( './testvtk/pointdata.txt' );
  tau = load ( './testvtk/celldata.txt' );

 %output_filename = './testvtk/verifycase/3.vtk';

  title = 'UVP for 1D advection-diffusion problem';
%
%  Have the data written to a file.
%
  twod_to_vtk ( nodes, elements, uvp,tau, output_filename, title );
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
