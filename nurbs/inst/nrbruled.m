function srf = nrbruled (crv1, crv2)

% NRBRULED: Construct a ruled surface between two NURBS curves, or a ruled volume between two NURBS surfaces.
% 
% Calling Sequence:
% 
%   srf = nrbruled(crv1, crv2)
% 
% INPUT:
% 
%   crv1	: First NURBS curve (or surface), see nrbmak.
% 
%   crv2	: Second NURBS curve (or surface), see nrbmak.
%
% OUTPUT:
% 
%   srf		: Ruled NURBS surface (or volume).
% 
% Description:
% 
%   Constructs a ruled surface between two NURBS curves. The ruled surface is
%   ruled along the V (or W) direction.
% 
% Examples:
% 
%   Construct a ruled surface between a semicircle and a straight line.
% 
%   cir = nrbcirc(1,[0 0 0],0,pi);
%   line = nrbline([-1 0.5 1],[1 0.5 1]);
%   srf = nrbruled(cir,line);
%   nrbplot(srf,[20 20]);
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2018 Rafael Vazquez
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

ndim = numel (crv1.order);
if (ndim ~= numel (crv2.order))
  error ('Both NURBS must be of the same kind (both curves or both surfaces)')
elseif (ndim == 3 || numel (crv2.order) == 3)
  error ('The input cannot be a NURBS volume')
end

% ensure both curves have a common degree
d = max ([crv1.order, crv2.order]);
crv1 = nrbdegelev (crv1, d - crv1.order);
crv2 = nrbdegelev (crv2, d - crv2.order);

knt1 = crv1.knots;
knt2 = crv2.knots;
if (~iscell (knt1))
  knt1 = {knt1};
end
if (~iscell (knt2))
  knt2 = {knt2};
end

% merge the knot vectors, to obtain a common knot vector
ka = cell (1, ndim); kb = cell (1, ndim);
for idim = 1:ndim
  k1 = knt1{idim};
  k2 = knt2{idim};
  ku = unique ([k1 k2]);
  n = length (ku);
  ka{idim} = [];
  kb{idim} = [];
  for i = 1:n
    i1 = length (find (k1 == ku(i)));
    i2 = length (find (k2 == ku(i)));
    m = max (i1, i2);
    ka{idim} = [ka{idim} ku(i)*ones(1,m-i1)];
    kb{idim} = [kb{idim} ku(i)*ones(1,m-i2)];
  end
end

if (ndim == 1)
  crv1 = nrbkntins (crv1, ka{1});
  crv2 = nrbkntins (crv2, kb{1});
  knots = {crv1.knots, [0 0 1 1]};
  coefs(:,:,1) = crv1.coefs;
  coefs(:,:,2) = crv2.coefs;
else
  crv1 = nrbkntins (crv1, ka);
  crv2 = nrbkntins (crv2, kb);
  knots = {crv1.knots{:}, [0 0 1 1]};
  coefs(:,:,:,1) = crv1.coefs;
  coefs(:,:,:,2) = crv2.coefs;
end
srf = nrbmak (coefs, knots);

end

%!demo
%! pnts = [0.5 1.5 4.5 3.0 7.5 6.0 8.5;
%!         3.0 5.5 5.5 1.5 1.5 4.0 4.5;
%!         0.0 0.0 0.0 0.0 0.0 0.0 0.0];
%! crv1 = nrbmak (pnts,[0 0 0 1/4 1/2 3/4 3/4 1 1 1]);
%! crv2 = nrbtform (nrbcirc (4,[4.5;0],pi,0.0),vectrans([0.0 4.0 -4.0]));
%! srf = nrbruled (crv1,crv2);
%! nrbplot (srf,[40 20]);
%! title ('Ruled surface construction from two NURBS curves.');
%! hold off

%!demo
%! srf1 = nrbtestsrf;
%! srf2 = nrb4surf([0 0 -1], [10 0 -1], [0 10 -1], [10 10 -1]);
%! vol = nrbruled (srf1, srf2);
%! nrbkntplot (vol);
%! title ('Ruled volume construction from two NURBS surfaces')