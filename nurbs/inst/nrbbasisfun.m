function [B, id] = nrbbasisfun (points, nrb)

% NRBBASISFUN: Basis functions for NURBS
%
% Calling Sequence:
% 
%    B     = nrbbasisfun (u, crv)
%    B     = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun (pts, srf)
%   [B, N] = nrbbasisfun ({u, v, w}, vol)
%   [B, N] = nrbbasisfun (pts, vol)
%
%    INPUT:
%   
%      u   - parametric coordinates along u direction
%      v   - parametric coordinates along v direction
%      w   - parametric coordinates along w direction
%      pts - array of scattered points in parametric domain, array size: (ndim,num_points)
%      crv - NURBS curve
%      srf - NURBS surface
%      vol - NURBS volume
%   
%    If the parametric coordinates are given in a cell-array, the values
%     are computed in a tensor product set of points
%
%    OUTPUT:
%   
%      B - Value of the basis functions at the points
%          size(B)=[npts, prod(nrb.order)]
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   
%
%    Copyright (C) 2009 Carlo de Falco
%    Copyright (C) 2015 Jacopo Corno
%    Copyright (C) 2016 Rafael Vazquez
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

  if (   (nargin<2) ...
      || (nargout>2) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=numel(nrb.number))) ...
      || iscell(points) && numel(points) ~= numel(nrb.number) ...
      )
    error('Incorrect input arguments in nrbbasisfun');
  end
  
  if (~iscell (nrb.knots)) %% NURBS curve
    knt = {nrb.knots};
  else                     %% NURBS surface or volume
    knt = nrb.knots;
  end
  
  ndim = numel (nrb.number);
  w = reshape (nrb.coefs(4,:), [nrb.number 1]);
  
  for idim = 1:ndim
    if (iscell (points))
      pts_dim = points{idim};
    else
      pts_dim = points(idim,:);
    end
    sp{idim} = findspan (nrb.number(idim)-1, nrb.order(idim)-1, pts_dim, knt{idim});
    N{idim} = basisfun(sp{idim}, pts_dim, nrb.order(idim)-1, knt{idim});
    num{idim} = numbasisfun (sp{idim}, pts_dim, nrb.order(idim)-1, knt{idim}) + 1;
  end
  
  if (ndim == 1)
     id = num{1};
     B = reshape (w(num{1}), size(N{1})) .* N{1};
     B = bsxfun (@(x,y) x./y, B, sum(B,2));
%      B = B ./ sum(B,2);
  else
    if (iscell (points))
      npts_dim = cellfun (@numel, points);
      cumnpts = cumprod([1 npts_dim]);
      npts = prod (npts_dim);
      val_aux = 1;
      numaux = 1;
      cumorder = cumprod ([1 nrb.order]);
      cumnumber = cumprod ([1 nrb.number]);
      for idim = 1:ndim
        val_aux = kron (N{idim}, val_aux);
        num_dim = reshape (num{idim}, 1, npts_dim(idim), 1, nrb.order(idim));
        num_dim = repmat (num_dim, cumnpts(idim), 1, cumorder(idim), 1);
        
        num_prev = reshape (numaux, cumnpts(idim), 1, cumorder(idim), 1);
        num_prev = repmat (num_prev, 1, npts_dim(idim), 1, nrb.order(idim));
        numaux = sub2ind ([cumnumber(idim) nrb.number(idim)], num_prev, num_dim);
        numaux = reshape (numaux, cumnpts(idim+1), cumorder(idim+1));
      end
      B = reshape (val_aux, npts, prod (nrb.order));
      id = reshape (numaux, npts, prod (nrb.order));
      W = w(id);
      B = bsxfun (@(x,y) x./y, W.*B, sum (W .* B, 2));

    else
      npts = numel (points(1,:));
      B = zeros (npts, prod(nrb.order));
      id = zeros (npts, prod(nrb.order));
      
      local_num = cell (ndim, 1);
      for ipt = 1:npts
        val_aux = 1;
        for idim = 1:ndim
          val_aux = reshape (val_aux.' * N{idim}(ipt,:), 1, []);
%           val_aux2 = kron (N{idim}(ipt,:), val_aux);
          local_num{idim} = num{idim}(ipt,:);
        end
        [local_num{:}] = ndgrid (local_num{:});
        id(ipt,:) = reshape (sub2ind (nrb.number, local_num{:}), 1, size(id, 2));
        W = reshape (w(id(ipt,:)), size(val_aux));
        val_aux = W .* val_aux;
        B(ipt,:) = val_aux / sum (val_aux);
      end
    end
  end
end  

%!demo
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = [1 1 1 1];
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! B = nrbbasisfun (u, nrb);
%! xplot = sum(bsxfun(@(x,y) x.*y, B, x),2);
%! plot(xplot, B)
%! title('Cubic Bernstein polynomials')
%! hold off

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! B = nrbbasisfun (u, nrb);
%! xplot = sum(bsxfun(@(x,y) x.*y, B, x),2);
%!
%! yy = y; yy(1) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U); 
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,1), aux(1,:).', w(1)*aux(2,:).')
%! assert(B(:,1), w(1)*aux(2,:).', 1e-6)
%! 
%! yy = y; yy(2) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2, u);
%! %figure, plot(xplot, B(:,2), aux(1,:).', w(2)*aux(2,:).')
%! assert(B(:,2), w(2)*aux(2,:).', 1e-6)
%!
%! yy = y; yy(3) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,3), aux(1,:).', w(3)*aux(2,:).')
%! assert(B(:,3), w(3)*aux(2,:).', 1e-6)
%!
%! yy = y; yy(4) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! aux = nrbeval(nrb2,u);
%! %figure, plot(xplot, B(:,4), aux(1,:).', w(4)*aux(2,:).')
%! assert(B(:,4), w(4)*aux(2,:).', 1e-6)

%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! aux1 = linspace(0,1,m); aux2 = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1)});
%! u = rand (1, 30); v = rand (1, 10);
%! u = u - min (u); u = u / max (u);
%! v = v - min (v); v = v / max (v);
%! [B, N] = nrbbasisfun ({u, v}, nrb);
%! assert (sum(B, 2), ones(300, 1), 1e-6)
%! assert (all (all (B<=1)), true)
%! assert (all (all (B>=0)), true)
%! assert (all (all (N>0)), true)
%! assert (all (all (N <= prod (nrb.number))), true)
%! assert (max (max (N)),prod (nrb.number))
%! assert (min (min (N)),1)

%!test
%! p1 = 2;   p2 = 3;   p3 = 2;
%! n1 = 4;   n2 = 5;   n3 = 4;
%! Lx = 1;  Ly = 1;   Lz = 1; 
%! crv = nrbline([1 0], [2 0]);
%! nrb = nrbtransp (nrbrevolve (crv, [], [0 0 1], pi/2));
%! nrb = nrbextrude (nrb, [0 0 1]);
%! nrb = nrbdegelev (nrb, [p1-1, p2-2, p3-1]);
%! aux1 = linspace(0,1,n1); aux2 = linspace(0,1,n2); aux3 = linspace(0,1,n3);
%! nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1), aux3(2:end-1)});
%! 
%! u = rand (1, 12); v = rand (1, 10); w = rand (1, 15);
%! u = u - min (u); u = u / max (u);
%! v = v - min (v); v = v / max (v);
%! w = w - min (w); w = w / max (w);
%! [B, N] = nrbbasisfun ({u, v, w}, nrb);
%! assert (all(sum(B, 2) - ones(numel(u)*numel(v)*numel(w),1) < 1e-6))
%! assert (all (all (B <= 1)) == true)
%! assert (all (all (B >= 0)) == true)
%! assert (all (all (N >  0))  == true)
%! assert (all (all (N <= prod (nrb.number))) == true)
%! assert (max (max (N)) == prod (nrb.number))
%! assert (min (min (N))== 1)