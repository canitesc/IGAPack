function [ GeoMesh ] = init1DGeometryBeam
%creates the geometry parametrization for a linear beam
numPatches = 1;
GeoMesh = cell(numPatches, 1);
numberElements = 1;
p = 1;
L = 1; %length of the beam

%initialize the geometry as NURBS-toolbox object

%control points
coefs(1:2,1) = [0; 0];
coefs(1:2,2) = [L; 0];

%knot vector
knotU = [0 0 1 1];
%geo = nrbmak(coefs, knotU);
%nrbctrlplot(geo)
GeoMesh{1} = genGeoMesh1D( knotU, coefs, p, numberElements );





end

