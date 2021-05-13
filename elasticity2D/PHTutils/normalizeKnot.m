function [normKnot] = normalizeKnot(knotVector)
% Normalizes (linear scales) a knot vector so that the first entry is 0
% and the last entry is 1. 
% 
% Input
% -----
%   knotVector - (1D array) the input knot vector
%
% Output
% ------
%   normKnot - (1D array) the normalized knot vector
%
normKnot = 1/(knotVector(end)-knotVector(1))*(knotVector-knotVector(1));
end

