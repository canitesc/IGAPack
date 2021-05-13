function [tiling, tile, ver3d] = nrbspheretiling (topology, radius, center)
%
% NRBSPHERETILING:  Makes an array of NURBS patches representing a
%                   full or partial tiling of the sphere.
%                
%
% Calling Sequences:
%
%   [tiling, tile, ver3d] = nrbspheretiling
%   [tiling, tile, ver3d] = nrbspheretiling(topology)
%   [tiling, tile, ver3d] = nrbspheretiling(topology, [radius], [center])
% 
% INPUT:
% 
%   topology: String specifying the desired topology for the tiling.
%           Options are:    - 'cube' (default)
%                           - 'ico' for paired icosahedron (nonconforming)
%                           - 'rdode' for rhombic dodecahedron
%                           - 'rtria' for rhombic triacontahedron
%                           - 'dico' for deltoidal icositetrahedron
%                           - 'dhexe' for deltoidal hexecontahedron
%                           - 'octant' for a tiling of the first octant
%
%   radius: Radius of the sphere, default 1.0
%
%   center: Center of the sphere, default (0,0,0)
% 
% OUTPUT:
%
%   tiling: Structure array of NURBS objects representing the tiling of the
%           sphere.
%
%   tile:   NURBS object representing one tile of the unit sphere. The tile
%           will be a fourth-order rational Bezier patch.
% 
%   ver3d:  3-by-4 matrix with the coordinates of the 4 (ordered) vertices
%           on the unit sphere, ordered clockwise when viewed from outside the
%           sphere.
% 
% See also: nrbspheretile
%
% For more details, see:
%  Sander Dedoncker, Laurens Coox, Florian Maurin, Francesco Greco, Wim Desmet
%  Bézier tilings of the sphere and their applications in benchmarking multipatch isogeometric methods
%  Computer Meth. Appl. Mech. Engrg., 2018
%
% Copyright (C) 2017 Sander Dedoncker
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Complete input

if nargin < 1
    topology = 'cube';
end

if (nargin < 2 || isempty (radius))
  radius = 1;
elseif (radius <= 0)
  warning ('The radius should be positive. Taking radius equal to one')
  radius = 1;
end

if (nargin < 3 || isempty (center))
  center = [0 0 0];
end

%% Catch invalid input

if (nargin > 3)
    error('Too many input arguments');
end

%% Load tilings

switch lower (topology)
    case 'cube'
        [tile,ver3d] = CubeLoad;
        [tiling] = Multi_SphereFromCube(tile);

    case 'ico'
        [tile,ver3d] = IcoLoad;
        [tiling] = Multi_SphereFromIco(tile);
        disp ('Note that the meshes of the icosahedral tiling are non-conforming')

    case 'rdode'
        [tile,ver3d] = RDodeLoad;
        [tiling] = Multi_SphereFromRDode(tile);

    case 'rtria'
        [tile,ver3d] = RTriaLoad;
        [tiling] = Multi_SphereFromRTria(tile);

    case 'dico'
        [tile,ver3d] = DIcoLoad;
        [tiling] = Multi_SphereFromDIco(tile);

    case 'dhexe'
        [tile,ver3d] = DHexeLoad;
        [tiling] = Multi_SphereFromDHexe(tile);

    case 'octant'
        [tile,ver3d] = DIcoLoad;
        [tiling] = OctantFromDIco(tile);

    otherwise
        error ('Invalid topology')
end

for iptc = 1:numel(tiling)
  tiling(iptc) = nrbtform (nrbtform (tiling(iptc), vecscale(radius*[1 1 1])), vectrans(center));
end

end

%% Auxiliary functions

% Cube

function [tile,ver3d] = CubeLoad

ver3d = [1 -1 -1 1;1 1 -1 -1;-1 -1 -1 -1]/sqrt(3);
tile = nrbspheretile (ver3d);

end

function [cubesphere] = Multi_SphereFromCube(NB)

%multipatch structure
cubesphere(1) = NB;

%rotate the different patches
rot1 = vecroty(-pi/2);
cubesphere(2) = nrbtform(cubesphere(1),rot1);

rot2 = vecrotz(-pi/2);
for ii = 2:4
  cubesphere(ii+1) = nrbtform(cubesphere(ii),rot2);
end

rot3 = vecroty(pi);
cubesphere(6) = nrbtform(cubesphere(1),rot3);

end

% Paired Icosahedron

function [tile,ver3d] = IcoLoad

%% Define vertex points

gr = (1+sqrt(5))/2; % golden ratio

Q1 = [1/2;0;-gr/2]; % see en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
Q2 = [0;gr/2;-1/2];
Q3 = [-1/2;0;-gr/2];
Q4 = [0;-gr/2;-1/2];

R = norm(Q1); % sphere radius

ver3d = [Q1,Q2,Q3,Q4]/R; % rescale to unit sphere

%% Generate tile

tile = nrbspheretile(ver3d); % make tile

end

function [icosphere] = Multi_SphereFromIco(NB)

icosphere(1) = NB;

%rotate the different patches
rot1 = vecrot(2*pi/5,[0;(1+sqrt(5))/2/2;-1/2]);
rot2 = vecrotx(pi);
icosphere(6) = nrbtform(icosphere(1),rot2);
for ii = 1:4
  icosphere(ii+1) = nrbtform(icosphere(ii),rot1);
  icosphere(ii+6) = nrbtform(icosphere(ii+5),rot1);
end

end

% Rhombic dodecahedron

function [tile,ver3d] = RDodeLoad

%% Define vertex points

s2 = sqrt(2);
s3 = sqrt(3);

Q1 = [s2/2;0;-s2/2]; % see en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
Q2 = [0;s3/3;-s2/s3];
Q3 = [-s2/2;0;-s2/2];
Q4 = [0;-s3/3;-s2/s3];

R = norm(Q1); % sphere radius

ver3d = [Q1,Q2,Q3,Q4]/R; % rescale to unit sphere

%% Generate tile
       
tile = nrbspheretile(ver3d); % make tile

end

function [rdodesphere] = Multi_SphereFromRDode(NB)

%multipatch structure
rdodesphere(1) = NB;

%rotate the different patches
rot1 = vecrot(pi/2,[sqrt(2)/2;0;-sqrt(2)/2]);
rdodesphere(2) = nrbtform(rdodesphere(1),rot1);
rdodesphere(3) = nrbtform(rdodesphere(2),rot1);
rdodesphere(4) = nrbtform(rdodesphere(3),rot1);

rot2 = vecrot(pi/2,[-sqrt(2)/2;0;-sqrt(2)/2]);
rdodesphere(5) = nrbtform(rdodesphere(1),rot2);
rdodesphere(6) = nrbtform(rdodesphere(5),rot2);
rdodesphere(7) = nrbtform(rdodesphere(6),rot2);

rot3 = vecrotx(pi);
rdodesphere(8) = nrbtform(rdodesphere(1),rot3);
rdodesphere(9) = nrbtform(rdodesphere(2),rot3);
rdodesphere(10) = nrbtform(rdodesphere(4),rot3);
rdodesphere(11) = nrbtform(rdodesphere(5),rot3);
rdodesphere(12) = nrbtform(rdodesphere(7),rot3);

rdodesphere = rdodesphere([1 2 7 10 11 6 5 4 3 8 12 9]);

end

% Rhombic triacontahedron

function [tile,ver3d] = RTriaLoad

%% Define vertex points

gr = (1+sqrt(5))/2; % golden ratio

Q1 = [gr^2;0;-gr^3]; % see http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/RhombicTriaconta/RhombicTriaconta.pdf
Q2 = [0;gr;-gr^3];
Q3 = [-gr^2;0;-gr^3];
Q4 = [0;-gr;-gr^3];

R1 = norm(Q1); R2 = norm(Q2); R3 = norm(Q3); R4 = norm(Q4);% sphere radius

ver3d = [Q1/R1,Q2/R2,Q3/R3,Q4/R4]; % rescale to unit sphere

%% Generate tile

tile = nrbspheretile(ver3d); % make tile

end
function [rtriasphere] = Multi_SphereFromRTria(NB)

%multipatch structure
rtriasphere(1) = NB;

%rotate the different patches
gr = (1+sqrt(5))/2;

rot1 = vecrot(2*pi/5,[gr^2;0;-gr^3]);
for ii = 1:3
  rtriasphere(ii+1) = nrbtform(rtriasphere(ii),rot1);
end

rot2 = vecrot(2*pi/5,[-gr^2;0;-gr^3]);
for ii = 1:4
  rtriasphere(ii+4) = nrbtform(rtriasphere(ii),rot2);
end
for ii = 1:12
  rtriasphere(ii+8) = nrbtform(rtriasphere(ii+4),rot2);
end

rot3 = vecroty(pi);
for ii = 1:5
  rtriasphere(20+ii) = nrbtform(rtriasphere(4*(ii-1)+1),rot3);
  rtriasphere(25+ii) = nrbtform(rtriasphere(4*(ii-1)+2),rot3);
end

end

% Deltoidal icositetrahedron

function [tile,ver3d] = DIcoLoad

%% Define vertex points

s3 = sqrt(3);

Q1 = [0;0;-1]; % quadrisected cube tile
Q2 = [s3/3;0;-s3/3];
Q3 = [s3/3;s3/3;-s3/3];
Q4 = [0;s3/3;-s3/3];

R1 = norm(Q1); R2 = norm(Q2); R3 = norm(Q3); R4 = norm(Q4);% sphere radius

ver3d = [Q1/R1,Q2/R2,Q3/R3,Q4/R4]; % rescale to unit sphere

%% Generate tile

tile = nrbspheretile(ver3d); % make tile

end
function [dicosphere] = Multi_SphereFromDIco(NB)

%multipatch structure
dicosphere(1) = NB;

%rotate the different patches

rot1 = vecrotz(pi/2);
for ii = 1:3
  dicosphere(ii+1) = nrbtform(dicosphere(ii),rot1);
end

rot2 = vecrotx(pi/2);
for ii = 1:12
  dicosphere(ii+4) = nrbtform(dicosphere(ii),rot2);
end

rot3 = vecroty(pi/2);
rot4 = vecroty(-pi/2);
for ii = 1:4
  dicosphere(16+ii) = nrbtform(dicosphere(ii),rot3);
  dicosphere(20+ii) = nrbtform(dicosphere(ii),rot4);
end

end

% Deltoidal hexecontahedron

function [tile,ver3d] = DHexeLoad

%% Define vertex points

gr = (1+sqrt(5))/2; % golden ratio

% Icosahedron
A1 = [1;0;-gr]; % see en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
A2 = [0;gr;-1];
A3 = [-1;0;-gr];

% Trisection
Q1 = A1;
Q2 = (A1+A2)/2;
Q3 = (A1+A2+A3)/3;
Q4 = (A3+A1)/2;

R1 = norm(Q1); R2 = norm(Q2); R3 = norm(Q3); R4 = norm(Q4);% sphere radius

ver3d = [Q1/R1,Q2/R2,Q3/R3,Q4/R4]; % rescale to unit sphere

% Rotation
% rotm = axang2rotm([ver3d(:,1)',-pi/5]);
rotm = vecrot(-pi/5,ver3d(:,1)');
rotm = rotm(1:3,1:3);
ver3d = rotm*ver3d;

%% Generate tile

tile = nrbspheretile(ver3d); % make tile

end
function [dhexesphere] = Multi_SphereFromDHexe(NB)

%multipatch structure
dhexesphere(1) = NB;

%make a paired ico patch
gr = (1+sqrt(5))/2;
A1 = [1;0;-gr]; % see en.wikipedia.org/wiki/Regular_icosahedron#Cartesian_coordinates
A2 = [0;gr;-1];
A3 = [-1;0;-gr];
A4 = (A1+A2+A3)/3;

rot0 = vecrot(pi/5,A1); %rotate first tile
dhexesphere(1) = nrbtform(dhexesphere(1),rot0);

rot1 = vecrot(2*pi/3,A4);
dhexesphere(2) = nrbtform(dhexesphere(1),rot1);
dhexesphere(3) = nrbtform(dhexesphere(2),rot1);

rot2 = vecrot(-2*pi/5,A1);
for ii = 1:3
  dhexesphere(ii+3) = nrbtform(dhexesphere(ii),rot2);
end

%rotate the 10 different patches
rot3 = vecrot(2*pi/5,A2);
for ii = 1:24
  dhexesphere(ii+6) = nrbtform(dhexesphere(ii),rot3);
end

rot4 = vecrotx(pi);
for ii = 1:30
  dhexesphere(30+ii) = nrbtform(dhexesphere(ii),rot4);
end

end


% Octant
function [octant] = OctantFromDIco(NB)

%rotate the different patches

rot1 = vecroty(-pi/2);
octant(1) = nrbtform(NB,rot1);

rot2 = vecrotx(pi/2);
octant(2) = nrbtform(NB,rot2);

rot3 = vecrotz(-pi/2);
temp = nrbtform(NB,rot3);
temp = nrbtform(temp,rot2);
octant(3) = nrbtform(temp,rot2);

end

%!demo
%! tiling = nrbspheretiling('cube');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical cubic tiling')

%!demo
%! tiling = nrbspheretiling('ico');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical icosahedral tiling')

%!demo
%! tiling = nrbspheretiling('rdode');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical rhombic dodecahedral tiling')

%!demo
%! tiling = nrbspheretiling('rtria');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical rhombic triacontahedral tiling')

%!demo
%! tiling = nrbspheretiling('dico');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical deltoidal icositetrahedral tiling')

%!demo
%! tiling = nrbspheretiling('dico');
%! figure
%! hold on
%! for ii = 1:length(tiling)
%!  nrbkntplot(tiling(ii), 10)
%! end
%! title('Spherical deltoidal hexecontahedral tiling')
