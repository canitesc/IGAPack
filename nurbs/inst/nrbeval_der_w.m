function der = nrbeval_der_w (nrb, i, points)
%
% NRBEVAL_DER_W: Compute the derivatives of a NURBS object at the point u 
% with respect to the weight of the i-th control point.
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
%   der - Derivatives.
%          size(der) = [3, numel(u)] for curves
%          or          [3, numel(u)*numel(v)] for surfaces
%          or          [3, numel(u)*numel(v)*numel(w)] for volumes
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
  
  if (iscell(points))
    npts = prod (cellfun (@numel, points));
  else
    npts = size (points, 2);
  end
    
  der = zeros (3, npts);
  [evalu, den] = nrbeval (nrb, points);
  [N, I] = nrbbasisfun (points, nrb);
  if (iscell (points))
    evalu = reshape (evalu, size(evalu, 1), []);
    den = reshape (den, 1, []);
  end
%   if (numel (nrb.number) == 1) % 1D
%     I = I + 1; % id is 0-based
%   end
  
  [ii, jj, kk] = ind2sub (nrb.number, i);
  w_i = nrb.coefs(4,ii,jj,kk);
  P_i = nrb.coefs(1:3,ii,jj,kk) ./ w_i;
  
  for ipnt = 1:npts
    [is, loc] = ismember (i, I(ipnt,:));
    if (is)
      der(:,ipnt) = N(ipnt,loc) ./ w_i .* P_i - evalu(:,ipnt) .* N(ipnt,loc) ./ w_i ./ den(ipnt);
    end
  end
    
end

%!test % 1D
%! nrb = nrbkntins (nrbcirc (1, [0 0], 0, pi/2), .5);
%! u = linspace (0, 1, 11);
%! delta_w = .01;
%! n = nrb.number;
%! der_ex = zeros (3, numel (u), n);
%! der_fd = zeros (3, numel (u), n);
%! for iw = 1:n
%!   new_w1   = nrb.coefs (4, iw) + delta_w;
%!   new_w2   = nrb.coefs (4, iw) - delta_w;
%!   nrb1 = nrbmodw (nrb, new_w1, iw);
%!   nrb2 = nrbmodw (nrb, new_w2, iw);
%!   der_ex(:,:,iw) = nrbeval_der_w (nrb, iw, u);
%!   p2 = nrbeval (nrb2, u);
%!   p1 = nrbeval (nrb1, u);
%!   der_fd(:,:,iw) = -(p2 - p1) ./ (2*delta_w);  
%! end
%! error = max (abs (der_ex(:) - der_fd(:)));
%! assert (error < 1.e-4)
%!
%!test %2D
%! crv = nrbline([1 0], [2 0]);
%! nrb = nrbtransp (nrbrevolve (crv, [], [0 0 1], pi/2));
%! new_knots = linspace (1/9, 8/9, 8);
%! nrb = nrbkntins (nrb, {new_knots, new_knots});
%! u = linspace (0, 1, 5);
%! v = u;
%! delta_w = .01;
%! n = nrb.number(1) * nrb.number(2);
%! der_ex = zeros (3, numel(u)* numel(v), n);
%! der_fd = zeros (3, numel(u)* numel(v), n);
%! for iw = 1:prod(nrb.number)
%!   new_w1   = nrb.coefs (4, iw) + delta_w;
%!   new_w2   = nrb.coefs (4, iw) - delta_w;
%!   nrb1 = nrbmodw (nrb, new_w1, iw);
%!   nrb2 = nrbmodw (nrb, new_w2, iw);
%!   der_ex(:,:,iw) = nrbeval_der_w (nrb, iw, {u v});
%!   p2 = nrbeval (nrb2, {u v});
%!   p1 = nrbeval (nrb1, {u v});
%!   der_fd(:,:,iw) = reshape (-(p2 - p1) ./ (2*delta_w), 3, []);
%! end
%! error = max (abs (der_ex(:) - der_fd(:)));
%! assert (error < 1.e-5)
%! 
%!test % 3D
%! crv = nrbline([1 0], [2 0]);
%! nrb = nrbtransp (nrbrevolve (crv, [], [0 0 1], pi/2));
%! nrb = nrbextrude (nrb, [0 0 1]);
%! u = 0:.33:.99;
%! v = 0:.1:.9;
%! w = [.25 .5 .75];
%! delta_w = .01;
%! n = nrb.number(1) * nrb.number(2) * nrb.number(3);
%! der_ex = zeros (3, numel(u)*numel(v)*numel(w), n);
%! der_fd = zeros (3, numel(u)*numel(v)*numel(w), n);
%! for iw = 1:prod(nrb.number)
%!   new_w1   = nrb.coefs (4, iw) + delta_w;
%!   new_w2   = nrb.coefs (4, iw) - delta_w;
%!   nrb1 = nrbmodw (nrb, new_w1, iw);
%!   nrb2 = nrbmodw (nrb, new_w2, iw);
%!   der_ex(:,:,iw) = nrbeval_der_w (nrb, iw, {u v w});
%!   p2 = nrbeval (nrb2, {u v w});
%!   p1 = nrbeval (nrb1, {u v w});
%!   der_fd(:,:,iw) = reshape (-(p2 - p1) ./ (2*delta_w), 3, []);
%! end
%! error = max (max (squeeze (max (abs (der_ex - der_fd)))));
%! assert (error < 1.e-4)
