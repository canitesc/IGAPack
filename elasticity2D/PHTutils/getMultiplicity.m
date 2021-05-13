function [uniqueKnots, mult] = getMultiplicity(knotVector)
% computes the multiplicities of a knot vector 
% Inputs
% ------
%   knotVector - 1D array of sorted knot vectors
%
% Outputs
% -------
%   uniqueKnots - 1D array of unique knot vectors
%   mult - 1D array with the multiplicity of each knot vector in uniqueKnots

uniqueKnots = unique(knotVector);
mult = zeros(1,length(uniqueKnots));
for i=1:length(uniqueKnots)
    mult(i) = length(find(knotVector==uniqueKnots(i)));    
end
assert(length(knotVector)==sum(mult))


