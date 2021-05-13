function der = nrbeval_der_p (nrb, i, points)
%
% NRBEVAL_DER_P: Compute the derivative of a NURBS object at a given point 
% with respect to the coordinates of the i-th control point.
%
% Calling Sequence:
% 
%   der = nrbeval_der_p (crv, i, u);
%   der = nrbeval_der_p (srf, i, p);
%   der = nrbeval_der_p (srf, i, {u v});
%   der = nrbeval_der_p (vol, i, p);
%   der = nrbeval_der_p (vol, i, {u v w});
%
% INPUT:
% 
%   crv	   - NURBS curve.
%   srf	   - NURBS surface.
%   vol	   - NURBS volume.
%   i      - Index of the control point.
%   u or p(1,:,:)  - parametric points along u direction
%   v or p(2,:,:)  - parametric points along v direction
%   w or p(3,:,:)  - parametric points along w direction
%
% OUTPUT:
%
%   der - Derivative.
%          size(der) = numel(u) for curves
%          or          numel(u)*numel(v) for surfaces
%          or          numel(u)*numel(v)*numel(w) for volumes
%   
% Copyright (C) 2015 Jacopo Corno
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
%
  
  [N, id] = nrbbasisfun (points, nrb);
  
  der = zeros (1, size(N, 1));
  for k = 1:numel (der)
    [is, loc] = ismember (i, id(k,:)); % id is 1-based

    if (is)
      der(k) = N(k,loc);
    else
      der(k) = 0;
    end
  end

end

%!test %% 1D
%! nrb =   nrbkntins (nrbcirc (1, [0 0], 0, pi/2), .5);
%! u = 0:.1:.9;
%! index = 1:nrb.number;
%! e = zeros (numel (u), numel (index), 1);
%! for jj = 1:numel (index)
%!   deltap = .1 * rand (3, 1);
%!   nrb2 = nrbmodp (nrb, deltap, index(jj));
%!   der_ex = nrbeval_der_p (nrb, index(jj), u);
%!   p2 = nrbeval (nrb2, u);
%!   p1 = nrbeval (nrb, u);
%!   der_fd = (p2 - p1) ./ deltap;
%!   e(:,jj) = sqrt (sum ((repmat (der_ex, 3, 1) - der_fd).^2, 1));
%! end
%! assert (max(e(:)) < 1.e-8);
%! 
%!test %% 2D
%! crv = nrbline([1 0], [2 0]);
%! nrb = nrbtransp (nrbrevolve (crv, [], [0 0 1], pi/2));
%! new_knots = linspace (1/9, 8/9, 8);
%! nrb = nrbkntins (nrb, {new_knots, new_knots});
%! u = 0:.1:.9;
%! v = u;
%! e = zeros (nrb.number(1) * nrb.number(2), numel (u), numel (v));
%! for index = 1:prod(nrb.number)
%!   deltap = .1 * rand (3, 1);
%!   nrb2 = nrbmodp (nrb, deltap, index);
%!   der_ex = nrbeval_der_p (nrb, index, {u v});
%!   p2 = nrbeval (nrb2, {u v});
%!   p1 = nrbeval (nrb, {u v});
%!   der_fd = (p2 - p1) ./ deltap;
%!   der_ex = reshape (repmat (der_ex, 3, 1), size(der_fd));
%!   e(index,:,:) = sqrt (sum ((der_ex - der_fd).^2, 1));
%! end
%! assert (max(e(:)) < 1.e-8)
%! 
%!test %% 3D
%! crv = nrbline([1 0], [2 0]);
%! nrb = nrbtransp (nrbrevolve (crv, [], [0 0 1], pi/2));
%! nrb = nrbextrude (nrb, [0 0 1]);
%! u = 0:.1:.9;
%! v = u;
%! w = u;
%! e = zeros (nrb.number(1) * nrb.number(2) * nrb.number(3), numel(u), numel(v), numel(w));
%! for index = 1:prod(nrb.number)
%!   deltap = .1 * rand (3, 1);
%!   nrb2 = nrbmodp (nrb, deltap, index);
%!   der_ex = nrbeval_der_p (nrb, index, {u v w});
%!   p2 = nrbeval (nrb2, {u v w});
%!   p1 = nrbeval (nrb, {u v w});
%!   der_fd = (p2 - p1) ./ deltap;
%!   der_ex = reshape (repmat (der_ex, 3, 1), size (der_fd));
%!   e(index,:,:,:) = sqrt (sum ((der_ex - der_fd).^2, 1));
%! end
%! assert (max (e(:)) < 1.e-8);
