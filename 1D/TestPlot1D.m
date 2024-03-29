%script TestPlot1D.m
% plot of basis functions before and after cross insertion
%pre-processing

clear
close all

addpath('./utils1D')
addpath('./ExampleData')
addpath('../nurbs/inst')

GeoMesh = init1DGeometryBeam;


numberElements = 4; %initial number of non-empty knot-spans; this should be an even number!

target_rel_error = 1e-10;
targetScale = 0.5;

p = 3; %polynomial degree


numPatches = length(GeoMesh);
PHTelem = cell(1,numPatches);
[ PHTelem{1}, dimBasis ] = initPHTmesh1DNoRefine( p,numberElements);
plotBasisPHTMesh1D( PHTelem, GeoMesh, p )
[ PHTelem{1}, dimBasis ] = crossInsert1D( PHTelem{1}, 3, dimBasis, p );
plotBasisPHTMesh1D( PHTelem, GeoMesh, p )