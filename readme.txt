IGAPack - A toolbox for isogeometric analysis using high-order PHT-splines

Several benchmark and test problems in 2D and 3D elasticity are implemented.
The main files are located in the top level of the elasticity2D and 
elasticity 3D directories (e.g. EdgeCrackMultiPatch.m or Hemisphere.m). The
examples have been tested with Matlab 2015b or later. For 3D, the deCasteljau
algorithm is implemented via mex files obtained with the Matlab Coder toolbox.
A pure MATLAB function (deCasteljau3DTP.m) is also available, but it will
probably be slower. The 3D solutions are outputed in the VTU format which is
viewable with Paraview (www.paraview.org) or similar software.

This program is distributed with version 1.3.8 of the NURBS toolbox. See the
nurbs directory for detailed licensing and usage information. 