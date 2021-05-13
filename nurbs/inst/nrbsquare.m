function srf = nrbsquare (corner, lengthx, lengthy, varargin)
%
% NRBSQUARE: create the NURBS surface for a square.
%
% Calling Sequence:
% 
% srf = nrbsquare (corner, lengthx, lengthy);
% srf = nrbsquare (corner, lengthx, lengthy, degree);
% srf = nrbsquare (corner, lengthx, lengthy, degree, nsub);
%
% INPUT:
%   corner  : the coordinates of the bottom left corner of the square. 
%   lenghtx : the length along the x direction.
%   lenghty : the length along the y direction.
%   degree  : the degree of the NURBS surface, in each direction.
%   nsub    : the number of subdivision of the NURBS surface, in each direction.
%
% OUTPUT:
%   srf     : the NURBS surface.
%
%    Copyright (C) 2016 Jacopo Corno, Rafael Vazquez
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

if (isempty (corner))
  corner = [0 0];
end

nsub = [1 1];
degree = [1 1];

if (numel (varargin) >= 1)
  if (numel (varargin{1}) == 1)
    degree = [varargin{1} varargin{1}];
  elseif (numel (varargin{1}) == 2)
    degree = varargin{1};
  else
    error ('The degree vector should provide the degree in each direction (two values).');
  end

  if (numel (varargin) == 2)
    if (numel (varargin{2}) == 1)
      nsub = [varargin{2} varargin{2}];
    elseif (numel (varargin{2}) == 2)
      nsub = varargin{2};
    else
      error ('The nsub vector should provide the number of intervals in each direction (two values).');
    end
  end
end

srf = nrb4surf (corner, corner+[lengthx 0], corner+[0 lengthy], corner+[lengthx lengthy]);
srf = nrbdegelev (srf, degree-[1 1]);
[~,~,new_knots] = kntrefine (srf.knots, nsub-1, degree, degree-[1 1]);
srf = nrbkntins (srf, new_knots);

end

%!test
%! srf = nrbsquare ([], 1, 2, 2, 4);
%! assert (srf.order, [3 3]);
%! knt = [0 0 0 1/4 1/2 3/4 1 1 1];
%! assert (srf.knots, {knt knt})
%! x = linspace (0, 1, 100);
%! [X,Y] = ndgrid (x, x);
%! vals = nrbeval (srf, {x x});
%! assert (squeeze(vals(1,:,:)), X, 1e-15);
%! assert (squeeze(vals(2,:,:)), 2*Y, 1e-15);
