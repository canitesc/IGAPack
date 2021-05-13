function idx = nrbnumbasisfun (points, nrb)
%
% NRBNUMBASISFUN:  Numbering of basis functions for NURBS
%
% Calling Sequence:
% 
%   N      = nrbnumbasisfun (u, crv)
%   N      = nrbnumbasisfun ({u, v}, srf)
%   N      = nrbnumbasisfun (p, srf)
%   N      = nrbnumbasisfun ({u, v, w}, vol)
%   N      = nrbnumbasisfun (p, vol)
%
%    INPUT:
%   
%      u or p(1,:,:)  - parametric points along u direction
%      v or p(2,:,:)  - parametric points along v direction
%      w or p(3,:,:)  - parametric points along w direction
%      crv - NURBS curve
%      srf - NURBS surface
%      vol - NURBS volume
%   
%    OUTPUT:
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == [npts, prod(nrb.order)]
%
%    Copyright (C) 2009 Carlo de Falco
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
      || (nargout>1) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=numel(nrb.number))) ...
      )
    error('Incorrect input arguments in nrbnumbasisfun');
  end


  if (~iscell(nrb.knots))          %% NURBS curve
    
    iv  = findspan (nrb.number-1, nrb.order-1, points, nrb.knots);
    idx = numbasisfun (iv, points, nrb.order-1, nrb.knots);

  else
    ndim = numel (nrb.number);
    if (iscell (points))
      for idim = 1:ndim
        pts_dim = points{idim};
        sp{idim} = findspan (nrb.number(idim)-1, nrb.order(idim)-1, pts_dim, nrb.knots{idim});
%         N{idim} = basisfun(sp{idim}, pts_dim, nrb.order(idim)-1, nrb.knots{idim});
        num{idim} = numbasisfun (sp{idim}, pts_dim, nrb.order(idim)-1, nrb.knots{idim}) + 1;
      end
      npts_dim = cellfun (@numel, points);
      cumnpts = cumprod([1 npts_dim]);
      npts = prod (npts_dim);

      numaux = 1;
      cumorder = cumprod ([1 nrb.order]);
      cumnumber = cumprod ([1 nrb.number]);
      for idim = 1:ndim
        num_dim = reshape (num{idim}, 1, npts_dim(idim), 1, nrb.order(idim));
        num_dim = repmat (num_dim, cumnpts(idim), 1, cumorder(idim), 1);
        
        num_prev = reshape (numaux, cumnpts(idim), 1, cumorder(idim), 1);
        num_prev = repmat (num_prev, 1, npts_dim(idim), 1, nrb.order(idim));
        numaux = sub2ind ([cumnumber(idim) nrb.number(idim)], num_prev, num_dim);
        numaux = reshape (numaux, cumnpts(idim+1), cumorder(idim+1));
      end
      idx = reshape (numaux, npts, prod (nrb.order));
      
    else
      for idim = 1:ndim
        pts_dim = points(idim,:);
        sp{idim} = findspan (nrb.number(idim)-1, nrb.order(idim)-1, pts_dim, nrb.knots{idim});
%         N{idim} = basisfun(sp{idim}, pts_dim, nrb.order(idim)-1, nrb.knots{idim});
        num{idim} = numbasisfun (sp{idim}, pts_dim, nrb.order(idim)-1, nrb.knots{idim}) + 1;
      end
      npts = numel (points(1,:));
      idx = zeros (npts, prod(nrb.order));
      
      local_num = cell (ndim, 1);
      for ipt = 1:npts
        for idim = 1:ndim
          local_num{idim} = num{idim}(ipt,:);
        end
        [local_num{:}] = ndgrid (local_num{:});
        idx(ipt,:) = reshape (sub2ind (nrb.number, local_num{:}), 1, size(idx, 2));
      end
    end
  end
  
end


%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! ikx = linspace(0,1,m); iky = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {ikx(2:end-1), iky(2:end-1)});
%! nrb.coefs (4,:,:) = nrb.coefs (4,:,:) + rand (size (nrb.coefs (4,:,:)));
%! u = rand (1, 30); v = rand (1, 10);
%! u = (u-min (u))/max (u-min (u));
%! v = (v-min (v))/max (v-min (v));
%! N = nrbnumbasisfun ({u, v}, nrb);
%! assert (all (all (N>0)), true)
%! assert (all (all (N <= prod (nrb.number))), true)
%! assert (max (max (N)), prod (nrb.number))
%! assert (min (min (N)), 1)