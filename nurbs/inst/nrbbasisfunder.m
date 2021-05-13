function varargout = nrbbasisfunder (points, nrb)

% NRBBASISFUNDER:  NURBS basis functions derivatives
%
% Calling Sequence:
% 
%   Bu          = nrbbasisfunder (u, crv)
%   [Bu, N]     = nrbbasisfunder (u, crv)
%   [Bu, Bv]    = nrbbasisfunder ({u, v}, srf)
%   [Bu, Bv, N] = nrbbasisfunder ({u, v}, srf)
%   [Bu, Bv, N] = nrbbasisfunder (pts, srf)
%   [Bu, Bv, Bw, N] = nrbbasisfunder ({u, v, w}, vol)
%   [Bu, Bv, Bw, N] = nrbbasisfunder (pts, vol)
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
%      Bu - Basis functions derivatives WRT direction u
%           size(Bu)=[npts, prod(nrb.order)]
%
%      Bv - Basis functions derivatives WRT direction v
%           size(Bv) == size(Bu)
%
%      Bw - Basis functions derivatives WRT direction w
%           size(Bw) == size(Bu)
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(Bu)
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
      || (nargout>4) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=numel(nrb.number))) ...
      || (~iscell(nrb.knots) && (nargout>2)) ...
      || iscell(points) && numel(points) ~= numel(nrb.number) ...
      )
    error('Incorrect input arguments in nrbbasisfunder');
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
    Nprime = basisfunder (sp{idim}, nrb.order(idim)-1, pts_dim, knt{idim}, 1);
    N{idim} = reshape (Nprime(:,1,:), numel(pts_dim), nrb.order(idim));
    Nder{idim} = reshape (Nprime(:,2,:), numel(pts_dim), nrb.order(idim));
    num{idim} = numbasisfun (sp{idim}, pts_dim, nrb.order(idim)-1, knt{idim}) + 1;
  end

  if (ndim == 1)
    B1 = reshape (w(num{1}), size(N{1})) .* N{1};
    W = sum (B1, 2);
    B2 = reshape (w(num{1}), size(N{1})) .* Nder{1};
    Wder = sum (B2, 2);

    B2 = bsxfun (@(x,y) x./y, B2, W);
    B1 = bsxfun (@(x,y) x.*y, B1, Wder./W.^2);
    B = B2 - B1;
    varargout{1} = B;
    varargout{2} = num{1};
  else
    id = nrbnumbasisfun (points, nrb);
    if (iscell (points))
      npts_dim = cellfun (@numel, points);
      npts = prod (npts_dim);
      val_aux = 1;
      val_ders = repmat ({1}, ndim, 1);
      for idim = 1:ndim
        val_aux = kron (N{idim}, val_aux);        
        for jdim = 1:ndim
          if (idim == jdim)
            val_ders{idim} = kron(Nder{jdim}, val_ders{idim});
          else
            val_ders{idim} = kron(N{jdim}, val_ders{idim});
          end
        end
      end
      
      B1 = w(id) .* reshape (val_aux, npts, prod(nrb.order));
      W = sum (B1, 2);
      for idim = 1:ndim
        B2 = w(id) .* reshape (val_ders{idim}, npts, prod(nrb.order));
        Wder = sum (B2, 2);        
        varargout{idim} = bsxfun (@(x,y) x./y, B2, W) - bsxfun (@(x,y) x.*y, B1, Wder ./ W.^2);
      end
    else
      npts = numel (points(1,:));
      B = zeros (npts, prod(nrb.order));
      Bder = repmat ({B}, ndim, 1);

      for ipt = 1:npts
        val_aux = 1;
        val_ders = repmat ({1}, ndim, 1);
        for idim = 1:ndim
          val_aux = reshape (val_aux.' * N{idim}(ipt,:), 1, []);
%           val_aux = kron (N{idim}(ipt,:), val_aux);
          for jdim = 1:ndim
            if (idim == jdim)
              val_ders{idim} = reshape (val_ders{idim}.' * Nder{jdim}(ipt,:), 1, []);
            else
              val_ders{idim} = reshape (val_ders{idim}.' * N{jdim}(ipt,:), 1, []);
            end
          end
        end
        wval = reshape (w(id(ipt,:)), size(val_aux));
        val_aux = val_aux .* wval;
        W = sum (val_aux);
        for idim = 1:ndim
          val_ders{idim} = val_ders{idim} .* wval;
          Wder = sum (val_ders{idim});
          Bder{idim}(ipt,:) = bsxfun (@(x,y) x./y, val_ders{idim}, W) - bsxfun (@(x,y) x.*y, val_aux, Wder ./ W.^2);
        end
      end
      varargout(1:ndim) = Bder(1:ndim);
    end
    if (nargout > ndim)
      varargout{ndim+1} = id;
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
%! [Bu, id] = nrbbasisfunder (u, nrb);
%! plot(u, Bu)
%! title('Derivatives of the cubic Bernstein polynomials')
%! hold off

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! [Bu, id] = nrbbasisfunder (u, nrb);
%! #plot(u, Bu)
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10), 

%!test
%! U = [0 0 0 0 1/2 1 1 1 1];
%! x = [0 1/4 1/2 3/4 1] ;
%! y = [0 0 0 0 0];
%! w = rand(1,5);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 300); 
%! [Bu, id] = nrbbasisfunder (u, nrb); 
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10)

%!test
%! p = 2;   q = 3;   m = 4; n = 5;
%! Lx  = 1; Ly  = 1; 
%! nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
%! nrb = nrbdegelev (nrb, [p-1, q-1]);
%! aux1 = linspace(0,1,m); aux2 = linspace(0,1,n);
%! nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1)});
%! nrb.coefs (4,:,:) = nrb.coefs(4,:,:) + rand (size (nrb.coefs (4,:,:)));
%! [Bu, Bv, N] = nrbbasisfunder ({rand(1, 20), rand(1, 20)}, nrb);
%! #plot3(squeeze(u(1,:,:)), squeeze(u(2,:,:)), reshape(Bu(:,10), 20, 20),'o')
%! assert (sum (Bu, 2), zeros(20^2, 1), 1e-10)

