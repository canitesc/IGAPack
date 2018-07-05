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

A description of the solving and error estimation procedure can be found in: https://doi.org/10.1016/j.cma.2017.08.032

Brief description of the main files

Elasticity 2D folder:

CantileverBeam.m: cantilever beam model from doi: 10.1016/j.finel.2008.01.010
EdgeCrackMultiPatch.m: plate with an edge crack model using 4 patches
ElasticRectangleMP.m: patch test for plate in tension with linear displacement, multiple patches can be used, Dirichlet boundary conditions
ElasticRectangeMPIso.m: patch test for plate in tension with linear displacement using isoparametric formulation, multiple patches can be used, Dirichlet boundary conditions
ElasticRectangleNeuMP.m: patch test for plate in tension with linear displacement, multiple patches can be used, traction boundary conditions
LshapedBracket.m: L-shape bracket geometry with 18 patches
LshapePanel.m: Benchmark L-shape panel with exact solution modeled with 3 patches
LshapePanelIso.m: Benchmark L-shape panel with exact solution modeled with 3 patches, isoparametric formulation
PlateHole.m: Benchmark plate with a hole problem modeled with 2 patches
PlateHoleC1.m: Benchmark plate with a hole problem, modeled with a single patch with a singular parametrization
PlateHoleIso.m: Benchmark plate with a hole problem modeled with 2 patches, isoparametric formulation
Spanner.m: Spanner problem modeled with 6 patches
SpannerIso.m: Spanner problem modeled with 6 patches, isoparametric formulation
ThickCylinder.m: Benchmark pressurized 2D cylinder problem

Elasticity 3D folder:

Blade.m: Turbine blade model
CantileverBar.m: Cantilever bar problem with 10 patches
CircularPlate.m: 3D disk model with vibration modes
ConnectingRod.m: Model of connecting rod
CubeWithHole.m: Benchmark cube with hole model with 1 patch
CubeWithHoleC0MP.m: Benchmark cube with hole model with 4 patches
EdgeCrackM2MP.m: Solid with an edge crack, mode II loading with exact solution
EdgeCrackM3MP.m: Solid with an edge crack, mode III loading
EdgeCrackMP.m: Solid with an edge crack, mode I loading with exact solution
ElasticBeamNeu.m: Cantilever bar with traction boundary conditions
ElasticCrackM2MP.m: Solid with an edge crack, mode II loading, simple boundary conditions
ElasticCrackM3MP.m: Solid with an edge crack, mode III loading, simple boundary conditions
ElasticCrackMP: Solid with an edge crack, mode I loading, simple boundary conditions
ElasticCubeMP: Patch test of solid with linear displacement, Dirichlet boundary conditions, multiple patches can be used
Hemisphere.m: Benchmark problem of 1/8th of a hollow sphere
HoreShoePHT.m: Horse shoe model with one patch
LshapedBracket3D.m: 3D L-Shaped bracket with 18 patches
PennyCrackMP.m: Benchmark problem of penny crack model
TbeamMP.m: Model of cantilever beam with varying thickness
