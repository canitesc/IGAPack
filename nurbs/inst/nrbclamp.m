function ccrv = nrbclamp (crv, k, xdim)

% NRBCLAMP: Compute the knot vector and control points of the clamped curve/surface.
%
% Calling Sequence:
% 
%   ccrv = nrbrclamp (crv)
%   ccrv = nrbrclamp (crv, k)
%   ccrv = nrbrclamp (crv, k, dim)
% 
% INPUT:
% 
%   crv	: unclamped NURBS curve or surface, see nrbmak.
%   k   : continuity desired afterclamping (from -1 up to p-1, -1 by default)
%   dim : dimension in which to clamp (all by default).
%
% OUTPUT:
% 
%   ccrv: NURBS curve with clamped knot vector, see nrbmak
% 
% Description:
% 
%   Clamps a curve or surface, using an open knot vector. Computes the new
%     knot vector and control points by knot insertion.
% 
%    Copyright (C) 2016 Monica Montardini, Filippo Remonato, Rafael Vazquez
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

if (iscell (crv.knots))
  knt = crv.knots;
  curve = false;
else
  knt = {crv.knots};
  curve = true;
end

ndim = numel (knt);
if (nargin < 2 || isempty(k))
  k = (-1) * ones (1, ndim);
end
if (nargin < 3)
  xdim = 1:ndim;
end

%if (iscell (crv.knots))
  if (numel(k) ~= ndim)
    k = k * ones(1, ndim);
  end

  new_knots = cell (1, ndim);
  for idim = xdim
    p  = crv.order(idim) - 1;
    U  = knt{idim};
    kk = k(idim);

    if (kk >= p)
      warning ('Taking the maximum k allowed, degree - 1')
      kk = p - 1;
    end

    n_ins_start(idim) = max (0, p - sum(U==U(p+1)) - kk);
    n_ins_end(idim) =  max (0, p - sum(U==U(end-p)) - kk);
    
    new_knots{idim} = [U(p+1)*ones(1,n_ins_start(idim)), U(end-p)*ones(1,n_ins_end(idim))];
  end

% Clamp, and remove unused coefficients and knots
  if (curve)
    ccrv = nrbkntins (crv, new_knots{1});
    ccrv.coefs = ccrv.coefs(:, n_ins_start+1 : end - n_ins_end);
    ccrv.knots = ccrv.knots(n_ins_start+1 : end - n_ins_end);
  else
    ccrv = nrbkntins (crv, new_knots);
    for idim = 1:ndim
      ccrv.knots{idim} = ccrv.knots{idim}(n_ins_start(idim)+1 : end - n_ins_end(idim));
      indices{idim} = n_ins_start(idim)+1 : ccrv.number(idim)-n_ins_end(idim);
    end
    ccrv.coefs = ccrv.coefs(1:4,indices{:});
  end
  
  ccrv.number = ccrv.number - n_ins_start - n_ins_end;
  
end

%!test
%! crv = nrbdegelev (nrbcirc (1, [], 0, pi/2), 2);
%! crv = nrbunclamp (crv, 3);
%! xx = linspace (0, 1, 20);
%! crv1 = nrbclamp (crv);
%! assert (crv1.knots, [0 0 0 0 0 1 1 1 1 1])
%! assert (nrbeval(crv, xx), nrbeval(crv1, xx), 1e-14)
%! crv1 = nrbclamp (crv, 2);
%! assert (crv1.knots, [-3 -2 -1 0 0 1 1 2 3 4])
%! assert (nrbeval(crv, xx), nrbeval(crv1, xx), 1e-14)

%!test
%! crv1 = nrbcirc(1,[],0,pi/4);
%! crv2 = nrbcirc(2,[],0,pi/4);
%! srf = nrbkntins (nrbdegelev (nrbruled(crv1, crv2), [3 2]), {0.25 []});
%! srf = nrbunclamp (srf, [4 2]);
%! srf1 = nrbclamp (srf);
%! xx = linspace(0,1,20);
%! assert(srf1.knots, {[0 0 0 0 0 0 0.2500 1 1 1 1 1 1]  [0 0 0 0 1 1 1 1]})
%! assert (nrbeval(srf, {xx xx}), nrbeval(srf1, {xx xx}), 1e-14);
%! srf1 = nrbclamp (srf, [3 1]);
%! assert (srf1.knots, {[-2 -1.75 -1 -0.75 0 0 0.25 1 1 1.25 2 2.25 3], [-2 -1 0 0 1 1 2 3]})
%! assert (nrbeval(srf, {xx xx}), nrbeval(srf1, {xx xx}), 1e-14);
%! srf1 = nrbclamp (srf, [], 2);
%! assert(srf1.knots, {[-2.75 -2 -1.75 -1 -0.75 0 0.25 1 1.25 2 2.25 3 3.25] [0 0 0 0 1 1 1 1]})
%! assert (nrbeval(srf, {xx xx}), nrbeval(srf1, {xx xx}), 1e-14);
