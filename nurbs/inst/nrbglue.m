function glued = nrbglue(nrb1, nrb2, side1, side2)
% 
% NRBGLUE: Glues two NURBS patches together with C^0-continuity at the
% interface. 
% 
% Calling Sequence:
% 
%   glued = nrbglue(nrb1, nrb2, [side1, side2])
% 
% INPUT:
% 
%   nrb1	: first NURBS struct.
%   nrb2	: second NURBS struct.
%   side1	: index of the boundary side of nrb1 to be glued (optional).
%   side2	: index of the boundary side of nrb2 to be glued (optional).
%
% OUTPUT:
% 
%   glued : Glued NURBS struct.
%
% The resulting NURBS struct has the same orientation as nrb1.
% 
% Description:
% 
%   The orientation of the NURBS boundaries:
% 
%                 ^ V direction
%                 |
%                 |     
%                 |              
%                 |             
%                 /---------> U direction
%                /
%               /
%              v  W direction
%
%   The indices of the NURBS boundaries (see also nrbextract): 
% 
%   1 : U = 0
%   2 : U = 1
%   3 : V = 0
%   4 : V = 1
%   5 : W = 0
%   6 : W = 1
%
%    Copyright (C) 2019 Ondine Chanon, Rafael Vazquez
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

interfaces = nrbmultipatch ([nrb1, nrb2]);
if (isempty (interfaces))
    error ('The patches cannot be glued together: non-compatible sides. Check degree and knot vectors')
end

if iscell(nrb1.knots)
    ndim = length(nrb1.knots);
else
    ndim = 1;
end

if (nargin == 4)
    if (side1 > 2*ndim || side2 > 2*ndim || side1 < 1 || side2 < 1)
        error ('Wrong number of the sides');
    end
    sides1 = [interfaces.side1];
    interface_index = find (sides1 == side1);
    if (isempty(interface_index) || interfaces(interface_index).side2 ~= side2)
        error ('The input sides are not correct')
    end
else
    if (numel (interfaces) > 1)
        error ('The two patches can be glued in more than one way (a ring?). Specify the sides')
    end
    side1 = interfaces.side1;
    side2 = interfaces.side2;
end
dir1 = ceil(side1/2);
dir2 = ceil(side2/2);

if (nrb1.order(dir1) < nrb2.order(dir2))
    deg_raise = zeros(1, ndim);
    deg_raise(dir1) = nrb2.order(dir2) - nrb1.order(dir1);
    nrb1 = nrbdegelev(nrb1, deg_raise);
    warning ('The order of nrb1 has been raised along the gluing direction')
elseif (nrb1.order(dir1) > nrb2.order(dir2))
    deg_raise = zeros(1, ndim);
    deg_raise(dir2) = nrb1.order(dir1) - nrb2.order(dir2);
    nrb2 = nrbdegelev(nrb2, deg_raise);
    warning ('The order of nrb2 has been raised along the gluing direction')
end

% Decide which patch will play the role of the first one
if (mod(side1,2) == 0)
    id1 = 1; id2 = 2;
else
    id1 = 2; id2 = 1;
end
if (mod(side1,2) == mod(side2,2))
    nrb2 = nrbreverse (nrb2, dir2);
end

% Make the gluing direction the first one
if (ndim > 1)
    dir_change = [dir1, setdiff(1:ndim, dir1)];
    nrb1 = nrbpermute (nrb1, dir_change);
    nrb2 = nrbpermute (nrb2, [dir2, setdiff(1:ndim, dir2)]);
end
nrbs = [nrb1, nrb2];
interfaces = nrbmultipatch (nrbs);

if (ndim == 1)
    nrbs(1).knots = {nrbs(1).knots};
    nrbs(2).knots = {nrbs(2).knots};
end

% Recompute, to simplify things
if (numel(interfaces) > 1)
    side = mod(side1-1,2) + 1;
    sides1 = [interfaces.side1];
    interfaces = interfaces(sides1 == side);
end

if (ndim == 2 && interfaces.ornt == -1)
    nrbs(2) = nrbreverse (nrbs(2), 2);
elseif (ndim == 3)
    if (interfaces.flag == -1)
        nrbs(2) = nrbpermute (nrbs(2), [1 3 2]);
    end
    if (interfaces.ornt1 == -1)
         nrbs(2) = nrbreverse (nrbs(2), 2);
    end
    if (interfaces.ornt2 == -1)
        nrbs(2) = nrbreverse (nrbs(2), 3);
    end
end

knots = nrbs(id1).knots;
knots{1} = [nrbs(id1).knots{1}(1:end-1), ...
        nrbs(id2).knots{1}(nrbs(id2).order(1)+1:end)+nrbs(id1).knots{1}(end)];
coefs = cat(2, nrbs(id1).coefs, nrbs(id2).coefs(:,2:end,:,:));

if (ndim == 1)
    knots = knots{1};
end

glued = nrbmak(coefs,knots);

% Recover the parametric directions of nrb1
if (ndim > 1)
  if (dir_change(1) ~= 3)
    glued = nrbpermute (glued, dir_change);
  else
    glued = nrbpermute (glued, [2 3 1]);
  end
end

%!demo
%! nrb1 = nrbline([0,0],[1,1]);
%! nrb2 = nrbline([1,1],[3,4]);
%! glued = nrbglue(nrb1, nrb2);
%! nrbkntplot(glued)
%! title ('Gluing of two straight lines')

%!demo
%! nrb1 = nrbsquare([0,0],1,1,[1,1],[1,2]);
%! nrb2 = nrbsquare([1,0],3,1,[1,1],[2,2]);
%! glued = nrbglue(nrb1, nrb2);
%! nrbkntplot(glued)
%! title ('Gluing of two bilinear patches')

%!demo
%! nrb1 = nrbsquare([0,0],1,1,[1,1],[1,2]);
%! nrb2 = nrbsquare([1,0],3,1,[1,1],[2,2]);
%! nrb1 = nrbextrude(nrb1,[0,0,1]);
%! nrb2 = nrbextrude(nrb2,[0,0,1]);
%! glued = nrbglue(nrb1, nrb2);
%! nrbkntplot(glued)
%! title ('Gluing of two trilinear patches')

