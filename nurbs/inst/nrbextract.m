function crvs = nrbextract(srf, sides)

%
% NRBEXTRACT: construct NURBS curves by extracting the boundaries of a NURBS surface, or NURBS surfaces by extracting the boundary of a NURBS volume.
% It only works for geometries constructed with open knot vectors. For a NURBS curve, 
% it returns two structures with the the boundary knots and control points.
% 
% Calling Sequence:
% 
%   crvs = nrbextract(surf, [sides]);
% 
% INPUT:
% 
%   surf        : NURBS surface or volume, see nrbmak.
%   sides       : the list of boundary sides to be extracted
% 
% OUTPUT: 
% 
%   crvs        : array of NURBS curves or NURBS surfaces extracted.
% 
% Description:
% 
%  Constructs either an array of four NURBS curves, by extracting the boundaries
%  of a NURBS surface, or an array of six surfaces, by extracting the boundaries
%  of a NURBS volume. The new entities are ordered in the following way
%
%    1: U = 0
%    2: U = 1
%    3: V = 0
%    4: V = 1
%    5: W = 0 (only for volumes)
%    6: W = 1 (only for volumes)
%
%    Copyright (C) 2010,2014,2015 Rafael Vazquez
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

if (nargin < 2)
  if (~iscell (srf.knots))
    ndim = 1;
  else
    ndim = numel(srf.knots);
  end
  sides = 1:2*ndim;
end

if (~iscell (srf.knots))
  crvs(1).knots = srf.knots(1);
  crvs(1).coefs = srf.coefs(:,1);
  crvs(2).knots = srf.knots(end);
  crvs(2).coefs = srf.coefs(:,end);
  crvs = crvs(sides);
  return
end

for idim = 1:numel(srf.knots)
  ord = srf.order(idim);
  if (srf.knots{idim}(1) ~= srf.knots{idim}(ord) || ...
      srf.knots{idim}(end) ~= srf.knots{idim}(end-ord+1))
    error ('nrbextract: only working for open knot vectors')
  end
end

if (numel (srf.knots) == 2)
  for ind = 1:2
    ind2 = mod (ind, 2) + 1;    %ind2 = [2 1];
    bnd1 = (ind - 1) * 2 + 1;
    bnd2 = (ind - 1) * 2 + 2;
    if (ind == 1)
      coefs1 = squeeze (srf.coefs(:,1,:));
      coefs2 = squeeze (srf.coefs(:,end,:));
    elseif (ind == 2)
      coefs1 = squeeze (srf.coefs(:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,end));
    end
    crvs(bnd1) = nrbmak (coefs1, srf.knots{ind2});
    crvs(bnd2) = nrbmak (coefs2, srf.knots{ind2});
  end  
elseif (numel (srf.knots) == 3)
  for ind = 1:3
    inds = setdiff (1:3, ind);
    bnd1 = (ind - 1) * 2 + 1;
    bnd2 = (ind - 1) * 2 + 2;
    if (ind == 1)
      coefs1 = squeeze (srf.coefs(:,1,:,:));
      coefs2 = squeeze (srf.coefs(:,end,:,:));
    elseif (ind == 2)
      coefs1 = squeeze (srf.coefs(:,:,1,:));
      coefs2 = squeeze (srf.coefs(:,:,end,:));
    elseif (ind == 3)
      coefs1 = squeeze (srf.coefs(:,:,:,1));
      coefs2 = squeeze (srf.coefs(:,:,:,end));
    end
    crvs(bnd1) = nrbmak (coefs1, {srf.knots{inds(1)} srf.knots{inds(2)}});
    crvs(bnd2) = nrbmak (coefs2, {srf.knots{inds(1)} srf.knots{inds(2)}});
  end
else
  error ('The entity is not a surface nor a volume')
end

crvs = crvs(sides);

end
