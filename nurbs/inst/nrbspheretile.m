function [tile] = nrbspheretile(ver3d)
%
% NRBSPHERETILE:    Makes a quadrilateral NURBS tile of the unit sphere
%                   from four vertex points.
%
% Calling Sequence:
% 
%   [tile] = nrbspheretile(ver3d)
% 
% INPUT:
% 
%   ver3d:  3-by-4 matrix with the coordinates of 4 (ordered) vertices on
%           the sphere. Vertices should be ordered clockwise when viewed
%           from outside the sphere. The tile should be convex and should
%           contain the south pole.
% 
% OUTPUT:
%
%   tile:   NURBS object representing the tile of the sphere. The tile
%           will be a quartic rational Bezier patch.
%
% See also: nrbspheretiling
%
% This function is based on the paper:
% J.E. Cobb, Tiling the Sphere with Rational Bezier patches, 1988
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

%% Catch invalid input

if nargin > 1
    error('Too many input arguments');
end

if size(ver3d,2) ~= 4 || size(ver3d,1) ~= 3
    error('Give a matrix of four column vectors.');
end

% ver3dnorm = round(sum(ver3d.^2,1),14);
ver3dnorm = sqrt (sum (ver3d.^2, 1));
if (any (abs (ver3dnorm - 1) > 1e-14))
    error('Vertices must lie on the unit sphere.')
end

%% Initialize variables
    
planes = zeros(3,4);   % matrix of a and b coeffs
circles = zeros(3,4);  % matrix of [Cu;Cv;R]
hecon = zeros(4,8);    % homogeneous coordinates of the external control points
hcon2d = zeros(4,3,3); % homogeneous coordinates for control points of the sgp tile
NM = zeros(25,25);     % coeff matrix of the collocation system
hcon3d = zeros(4,5,5); % homogeneous coordinates for control points of the sphere tile

knots_deg2 = [0 0 0 1 1 1];
knots_deg4 = [0 0 0 0 0 1 1 1 1 1];

%% Solve

ver3d = [ver3d,ver3d(:,1),ver3d(:,2)]; % repeat first and second point
ver2d = istf(ver3d); % representations of the points in the uv plane

% Find planes through great circle (ax+by+z=0 or ax+by=0)
% planes(:,i) = [ai;bi;1] or [ai;bi;0]
% Planes are defined by their inward (wrt the tile) pointing normal vectors

npchk = zeros(1,4);
spchk = zeros(1,4);

for ii = 1:4
    
    nvec1 = cross((ver3d(:,ii+1)-ver3d(:,ii)),ver3d(:,ii)); % inward pointing vector normal to the current plane
    nvec1 = nvec1./norm(nvec1);
    
    npchk(ii) = [0 0  1]*nvec1; % keeps track of location of north pole wrt the planes
    spchk(ii) = [0 0 -1]*nvec1; % keeps track of location of south pole wrt the planes
    
    if nvec1(3,1) == 0              % plane is vertical (ax+by=0)
        planes(:,ii) = nvec1(1:3,1);
    else                            % plane is not vertical (ax+by+z=0)
        planes(:,ii) = nvec1(1:3,1)./nvec1(3,1);
    end
    
    nvec2 = cross((ver3d(:,ii+2)-ver3d(:,ii+1)),ver3d(:,ii+1)); % inward pointing vector normal to the subsequent plane
    nvec2 = nvec2./norm(nvec2);
    
    tvec = cross(ver3d(:,ii+1),nvec2); % vector tangential to the next plane and normal to the intersection line
    tvec = tvec./norm(tvec); 
    
    if nvec1'*tvec<0 || any(isnan([nvec1;tvec]))
        error('Vertices are invalid. Check that: no two vertices are the same, the internal angles are smaller than 180 degrees, the vertices are ordered correctly.');
    end
    
end

if any(spchk < zeros(1,4))
    error('Tile does not contain the south pole, which may lead to a noninjective mapping. Consider rotating the vertices.')
end

if all(npchk >= zeros(1,4))
    error('Tile contains the north pole. Consider rotating the vertices.')
end

% Find circle arcs representing edges in the (u,v) plane. 
% R = 2*sqrt(1+a^2+b^2) and C = (-2a,-2b)

circles(1,:) = -2*planes(1,:);
circles(2,:) = -2*planes(2,:);
circles(3,:) = 2*sqrt(1 + planes(1,:).^2 + planes(2,:).^2);

% Find external control points in the z=-1 plane
for ii=1:4
    P1 = ver2d(:,ii);
    P3 = ver2d(:,ii+1);
    Q = (P1+P3)/2;
    C = circles(1:2,ii);
    R = circles(3,ii);
    
    if planes(3,ii) == 0 %[-ver3d(2,i),ver3d(1,i)]*(ver3d(1:2,i) - ver3d(1:2,i+1)) == 0 % plane is vertical (ax+by=0)
        P2 = Q; 
        w2 = 1; 
    else
        P2 = C+(Q-C)*R^2/norm(Q-C)^2; % |CP2| = |CQ|/(cos(alpha)^2)
        w2 = norm(Q-C)/R; % w2 = cos(alpha)
    end
    
    hecon(:,(ii-1)*2+1:ii*2)=[P1,P2*w2;-1,-1*w2;1,w2];% in 4D homogenous coordinates! Needed to construct patch properly
end

% Construct control net in the z=-1 plane

hcon2d(:,1,1) = hecon(:,1);
hcon2d(:,2,1) = hecon(:,2);
hcon2d(:,3,1) = hecon(:,3);
hcon2d(:,3,2) = hecon(:,4);
hcon2d(:,3,3) = hecon(:,5);
hcon2d(:,2,3) = hecon(:,6);
hcon2d(:,1,3) = hecon(:,7);
hcon2d(:,1,2) = hecon(:,8);
hcon2d(:,2,2) = [0;0;-1;1]; % central point

% Make patch in the z=-1 plane

sgp = nrbmak(hcon2d,{knots_deg2, knots_deg2});

% Calculate points corresponding to collocation points

xi = 0:1/4:1; % collocation points in parametric space
et = xi;

[p,w] = nrbeval(sgp,{xi,et});
hcol2d = [p;w];

% Project onto the sphere using the s-mapping

hcol3d = pstf(hcol2d);
     
% Evaluate the basis functions for the projective 3 space in the collocation points

bfxi = basisfun (findspan(4,4,xi,knots_deg4),xi,4,knots_deg4);
bfet = basisfun (findspan(4,4,et,knots_deg4),et,4,knots_deg4);

for ii = 1:5
    for jj = 1:5
       temp = bfxi(ii,:)'*bfet(jj,:);
       NM(ii + (jj-1)*5,:) = temp(:); % equation for collocation point k
    end
end

% Solve for each dimension separately

for ii = 1:4
    rhs = hcol3d(ii,:,:);
    temp = NM\rhs(:);
    hcon3d(ii,:,:) = reshape(temp,[1,5,5]);
end

% Make tile

tile = nrbmak(hcon3d,{knots_deg4, knots_deg4});

end

%% Auxiliary functions

function [ uv ] = istf( xyz )
%
% ISTF: s-mapping of the unit sphere in xyz-space to the uv plane. See
% Cobb, 1988
%
% Calling Sequences:
% 
%   [ xyz ] = istf( uv )
% 
% INPUT:
% 
%   xyz    : 3-by-n-by-m array with points in the xyz-space
%
% OUTPUT:
%
%   uv     : 2-by-n-by-m array with points in the uv-plane
% 

x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);

dim = size(xyz);
dim(1) = 2;

u = 2*x./(-z+1);
v = 2*y./(-z+1);

uv = [u;v];

uv = reshape(uv,dim);

end

function [ xyzw ] = pstf( uvw )
%
% PSTF: s-mapping of the projective uv plane to the unit sphere in 
% projective xyz-space. See Cobb, 1988
%
% Calling Sequences:
% 
%   [ xyzw ] = pstf( uvw )
% 
% INPUT:
% 
%   uvw     : 3-by-n-by-m array with points in the projective uv-plane
%
% OUTPUT:
%
%   xyzw    : 4-by-n-by-m array with points in the projective xyz-space
% 

u = uvw(1,:);
v = uvw(2,:);
w = uvw(end,:);

dim = size(uvw);
dim(1) = 4;

x = 4*u.*w;
y = 4*v.*w;
z = u.*u+v.*v-4*w.*w;
w2 = u.*u+v.*v+4*w.*w;

xyzw = [x;y;z;w2];

xyzw = reshape(xyzw,dim);

end

%!demo
%! vertices = [1 -1 -1 1;1 1 -1 -1;-1 -1 -1 -1]/sqrt(3);
%! tile = nrbspheretile(vertices);
%! figure
%! nrbkntplot(tile)
%! title('Spherical cube tile')

%!demo
%! vertices = [0 1 sqrt(2)/2 0;0 0 sqrt(2)/2 1; -1 0 0 0] ;
%! tile = nrbspheretile(vertices);
%! figure
%! nrbkntplot(tile)
%! title('Single patch octant tile')


