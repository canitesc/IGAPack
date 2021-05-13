function srf = bspinterpsurf (X, Y, Z, p, method)
%
% BSPINTERPSURF: B-Spline surface interpolation.
%
% Calling Sequence:
% 
%   srf = bspinterpsurf (Q, p, method);
%   
%    INPUT:
%   
%      X, Y, Z - grid of points to be interpolated. (See ndgrid)
%      p       - degree of the interpolating curve ([degree_x, degree_y]).
%      method  - parametrization method. The available choices are:
%                'equally_spaced'
%                'chord_length' (default)
%   
%    OUTPUT:
%   
%      srf - the B-Spline surface.
%   
%    See The NURBS book pag. 376 for more information. As of now only the
%    chord length method is implemented.
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

  if (nargin < 5 || isempty (method))
    method = 'chord_length';
  end
  
  [n, m] = size (X);
  Q = zeros (3, n, m);
  Q(1,:,:) = X;
  Q(2,:,:) = Y;
  Q(3,:,:) = Z;
  
  if (strcmpi (method, 'equally_spaced'))
    u = linspace (0, 1, n);
    v = linspace (0, 1, m);
  elseif (strcmp (method, 'chord_length'))
    u = zeros (m, n);
    for ii = 1:m
      d = sum (sqrt (sum (diff (squeeze(Q(:,:,ii))')'.^2,1)));
      u(ii,2:n) = cumsum (sqrt (sum (diff(Q(:,:,ii), [], 2).^2, 1)))/d;
%      for jj = 2:n-1 
%        u(ii,jj) = u(ii,jj-1) + norm (Q(:,jj,ii) - Q(:,jj-1,ii)) / d;
%      end
      u(ii,end) = 1;
    end
    u = mean (u);

    v = zeros (n, m);
    for ii = 1:n
      d = sum (sqrt (sum (diff (squeeze(Q(:,ii,:))')'.^2,1)));
      v(ii,2:m) = cumsum (sqrt (sum (diff(Q(:,ii,:), [], 3).^2, 1)))/d;
%      for jj = 2:m-1 
%        v(ii,jj) = v(ii,jj-1) + norm (Q(:,ii,jj) - Q(:,ii,jj-1)) / d;
%      end
      v(ii,end) = 1;
    end
    v = mean (v);
  end
  % TODO: implement centripetal method
  
  % Compute knot vectors
  knts{1} = zeros (1, n+p(1)+1);
  for jj = 2:n-p(1)
    knts{1}(jj+p(1)) = 1/p(1) * sum (u(jj:jj+p(1)-1));
  end
  knts{1}(end-p(1):end) = ones (1, p(1)+1);
  
  knts{2} = zeros (1, m+p(2)+1);
  for jj = 2:m-p(2)
    knts{2}(jj+p(2)) = 1/p(2) * sum (v(jj:jj+p(2)-1));
  end
  knts{2}(end-p(2):end) = ones (1, p(2)+1);
  
  % Interpolation
  R = zeros (size (Q));
  P = zeros (4, n, m);
  for ii = 1:m
    A = zeros (n, n);
    A(1,1) = 1;
    A(n,n) = 1;
    for jj=2:n-1
      span = findspan (n, p(1), u(jj), knts{1});
      A(jj,span-p(1)+1:span+1) = basisfun (span, u(jj), p(1), knts{1});
    end
    R(1,:,ii) = A \ squeeze(Q(1,:,ii))';
    R(2,:,ii) = A \ squeeze(Q(2,:,ii))';
    R(3,:,ii) = A \ squeeze(Q(3,:,ii))';
  end
  
  for ii = 1:n
    A = zeros (m, m);
    A(1,1) = 1;
    A(m,m) = 1;
    for jj=2:m-1
      span = findspan (m, p(2), v(jj), knts{2});
      A(jj,span-p(2)+1:span+1) = basisfun (span, v(jj), p(2), knts{2});
    end
    P(1,ii,:) = A \ squeeze(R(1,ii,:));
    P(2,ii,:) = A \ squeeze(R(2,ii,:));
    P(3,ii,:) = A \ squeeze(R(3,ii,:));
  end
  P(4,:,:) = ones (n, m);
  
  % Create B-Spline interpolant
  srf = nrbmak (P, knts);
  
end

%!demo
%! x = linspace (-3, 3, 40);
%! y = linspace (-3, 3, 40);
%! [X, Y] = meshgrid (x, y);
%! Z = peaks (X, Y);
%! 
%! srf1 = bspinterpsurf (X, Y, Z, [2 2], 'equally_spaced');
%! srf2 = bspinterpsurf (X, Y, Z, [2 2], 'chord_length');
%! figure
%! nrbkntplot(srf1)
%! title ('Approximation of the peaks functions, with the equally spaced method')
%! figure
%! nrbkntplot(srf2)
%! title ('Approximation of the peaks functions, with the chord length method')
