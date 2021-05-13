function [geomData] = genQuadWInclusion(vertices, center, rad)
%  Generate the geometry of a quadrilateral with inclusion domain 
%  Input: vertices: 4x2 array of the quadrilateral vertices in counter-clockwise
%                     order ([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], 
%                             [x4, y4, z4]])
%            center: coordinates of the center of the inclusion
%            rad: radius of the inclusion
%  Output: geomData : cell array of 4 NURBS structures corresponding to the
%          four patches

geomData = cell(1,4);

% set knot vectors
knot_u = [0, 0, 0, 1, 1, 1];
knot_v = [0, 0, 1, 1];

% set the locations of the control points scaled ot the rad and translated
% to center
w = sqrt(2)/2;

cpt_sw = [-w, -w, 0]*rad+center;
cpt_s = ([0, -1/w, 0]*rad+center)*w;
cpt_se = [w, -w, 0]*rad+center;
cpt_e = ([1/w, 0, 0]*rad+center)*w;
cpt_ne = [w, w, 0]*rad+center;
cpt_n = ([0, 1/w, 0]*rad+center)*w;
cpt_nw = [-w, w, 0]*rad+center;
cpt_w = ([-1/w, 0, 0]*rad+center)*w;

% set the location of the edge midpoints
vert_s = (vertices(1,:) + vertices(2,:))/2;
vert_e = (vertices(2,:) + vertices(3,:))/2;
vert_n = (vertices(3,:) + vertices(4,:))/2;
vert_w = (vertices(4,:) + vertices(1,:))/2;

% set the control points for the four patches
cptsSouth = [vertices(1,:); cpt_sw; vert_s; cpt_s; vertices(2,:); cpt_se];
cptsEast = [vertices(2,:); cpt_se; vert_e; cpt_e; vertices(3,:); cpt_ne];
cptsNorth = [vertices(3,:); cpt_ne; vert_n; cpt_n; vertices(4,:); cpt_nw];
cptsWest = [vertices(4,:); cpt_nw; vert_w; cpt_w; vertices(1,:); cpt_sw];

% set the weights
wgts = [1; 1; 1; w; 1; 1];

coefsSouth = reshape([cptsSouth, wgts]', 4, 2, 3);
coefsSouth = permute(coefsSouth, [1, 3, 2]);
coefsEast = reshape([cptsEast, wgts]', 4, 2, 3);
coefsEast = permute(coefsEast, [1, 3, 2]);
coefsNorth= reshape([cptsNorth, wgts]', 4, 2, 3);
coefsNorth = permute(coefsNorth, [1, 3, 2]);
coefsWest = reshape([cptsWest, wgts]', 4, 2, 3);
coefsWest = permute(coefsWest, [1, 3, 2]);

geomData{1} = nrbmak(coefsSouth, {knot_u, knot_v});
geomData{2} = nrbmak(coefsEast, {knot_u, knot_v});
geomData{3} = nrbmak(coefsNorth, {knot_u, knot_v});
geomData{4} = nrbmak(coefsWest, {knot_u, knot_v});


end

