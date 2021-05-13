function [geomData] = genQuad(vertices)
%  Generate the geometry of a quadrilateral
%  Input: vertices: 4x2 array of the quadrilateral vertices in counter-clockwise
%                     order ([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], 
%                            [x4, y4, z4]])
%            
%  Output: geomData :  NURBS structures corresponding to the bilinear quad
%                       patch


% set knot vectors
knot_u = [0, 0, 1, 1];
knot_v = [0, 0, 1, 1];

coefs = ones(4,2,2);
coefs(1:3,1,1) = vertices(1,:);
coefs(1:3, 2,1) = vertices(2,:);
coefs(1:3,2,2) = vertices(3,:);
coefs(1:3,1,2) = vertices(4,:);

geomData = nrbmak(coefs, {knot_u, knot_v});

end

