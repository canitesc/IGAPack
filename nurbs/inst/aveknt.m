function pts = aveknt (varargin)

% AVEKNT:  compute the knot averages (Greville points) of a knot vector
%
% Calling Sequence:
% 
%   pts = aveknt (knt, p)
%   pts = aveknt (nrb)
%   
%    INPUT:
%   
%      knt - knot sequence
%      p   - spline order (degree + 1)
%      nrb - NURBS structure (see nrbmak)
%   
%    OUTPUT:
%   
%      pts - average knots. If the input is a NURBS, it gives a cell-array,
%        with the average knots in each direction
%   
% See also: 
%
% Copyright (C) 2016 Rafael Vazquez
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

if (nargin == 1)
  if (isfield (varargin{1}, 'form'))
    nrb = varargin{1};
    knt = nrb.knots;
    order = nrb.order;
  else
    error ('The input should be a NURBS structure, or a knot vector and the order. See the help for details')
  end
elseif (nargin == 2)
  knt = varargin{1};
  order = varargin{2};
else
  error ('The input should be a NURBS structure, or a knot vector and the order. See the help for details')
end

onedim = false;
if (~iscell (knt))
  knt = {knt};
  onedim = true;
end

ndim = numel (knt);
pts = cell (ndim, 1);
for idim = 1:ndim
  if (numel (knt{idim}) < order(idim)+1)
    error ('The knot vector must contain at least p+2 knots, with p the degree')
  end
    
  knt_aux = repmat (knt{idim}(2:end-1), 1, order(idim)-1);
  knt_aux = [knt_aux(:); zeros(order(idim)-1, 1)];
  knt_aux = reshape (knt_aux, [], order(idim)-1);
  pts{idim} = sum (knt_aux.', 1) / (order(idim)-1);
  pts{idim} = pts{idim}(1:end-order(idim)+1);
end

if (onedim)
  pts = pts{1};
end

end

%!test
%! knt = [0 0 0 0.5 1 1 1];
%! pts = aveknt (knt, 3);
%! assert (pts - [0 1/4 3/4 1] < 1e-14)
%!
%!test
%! knt = {[0 0 0 0.5 1 1 1] [0 0 0 0 1/3 2/3 1 1 1 1]};
%! pts = aveknt (knt, [3 4]);
%! assert (pts{1} - [0 1/4 3/4 1] < 1e-14);
%! assert (pts{2} - [0 1/9 1/3 2/3 8/9 1] < 1e-14);
%!
%!test
%! nrb = nrb4surf([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbkntins (nrbdegelev (nrb, [1 2]), {[1/2] [1/3 2/3]});
%! pts = aveknt (nrb);
%! assert (pts{1} - [0 1/4 3/4 1] < 1e-14);
%! assert (pts{2} - [0 1/9 1/3 2/3 8/9 1] < 1e-14);
