function [ taux, tauy ] = parabolicTraction(coords, ~, params)
% outputs a parabolic traction 
% Input
% -----
%    coords - (array of length 2): the x and y coordinate of the evaluation
%               point
%
%    params - structure with fields: P - applied pressure, D - width of the
%             beam

taux = 0;
tauy = -params.P/2*(params.D^2/4-coords(2)^2);

