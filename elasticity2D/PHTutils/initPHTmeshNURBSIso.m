function [ PHTelem, controlPts, dimBasis, quadList ] = initPHTmeshNURBSIso(nrbList)
% Initialize the PHT geometry on coarse mesh from a list of NURBS toolbox 
% patches. It is  assumed that all the patches have at least C^1 continuity
% (no C^0 lines). The mesh is degree elevated to p=q=3 and  all knot-spans
% are uniformly refined once to create the recovery patches which are 
% stored in quadList.
% 
%   Inputs
%   ------
%   nrbList : cell array of NURBS toolbox structures containing the geometry
%             information
%
%   Outputs
%   -------
%   PHTelem    : cell array of PHTelem strctures containing the mesh information
%   controlPts : cell array of control points in each patch 
%   dimBasis   : size of the basis in each patch
%   quadlist   : cell array of recovery patches


numPatches = length(nrbList);
PHTelem = cell(1, numPatches);
controlPts = cell(1, numPatches);
dimBasis = zeros(1, numPatches);
quadList = cell(1, numPatches);

for i=1:numPatches
    knotU = nrbList{i}.knots{1};
    knotV = nrbList{i}.knots{2};
    numberElementsU = length(unique(knotU))-1;
    numberElementsV = length(unique(knotV))-1;
    p = 3;
    q = 3;
    [controlPts{i}, PHTelem{i}, dimBasis(i), quadList{i}] = genControlPts( knotU,...
        knotV, nrbList{i}.coefs, p, q, numberElementsU, numberElementsV );
    
end
