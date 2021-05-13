function nrbkntplot (nurbs, nsub)

% NRBKNTPLOT: Plot a NURBS entity with the knots subdivision.
% 
% Calling Sequence:
% 
%   nrbkntplot(nurbs)
%   nrbkntplot(nurbs, npnts)
% 
% INPUT:
% 
%   nurbs: NURBS curve, surface or volume, see nrbmak.
%   npnts: Number of evaluation points, for a surface or volume, a row 
%       vector with the number of points along each direction.
%
% Example:
%
%   Plot the test surface with its knot vector
%
%   nrbkntplot(nrbtestsrf)
%
% See also:
% 
%   nrbctrlplot
%
%    Copyright (C) 2011, 2012, 2021 Rafael Vazquez
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

if (nargin < 1)
  error ('nrbkntplot: Need a NURBS to plot!');
end

% Default values
light='on';
cmap='summer';

colormap (cmap);

hold_flag = ishold;

if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
   if (nargin < 2)
     nsub = [50 50];
   elseif (numel(nsub) == 1)
     nsub = [nsub nsub];
   end
   nrbplot (nurbs, nsub, 'light', light, 'colormap', cmap);
   hold on

   % And plot the knots
   knt1 = unique (nurbs.knots{1}(nurbs.order(1):end-nurbs.order(1)+1));
   knt2 = unique (nurbs.knots{2}(nurbs.order(2):end-nurbs.order(2)+1));
   p1 = nrbeval (nurbs, {knt1, linspace(knt2(1),knt2(end),nsub(2)+1)});
   p2 = nrbeval (nurbs, {linspace(knt1(1),knt1(end),nsub(1)+1), knt2});

  if (any (nurbs.coefs(3,:)))
    % surface in a 3D space
    for ii = 1:numel(knt1)
      plot3 (squeeze(p1(1,ii,:)), squeeze(p1(2,ii,:)), squeeze(p1(3,ii,:)),'k');
    end
    for ii = 1:numel(knt2)
      plot3 (squeeze(p2(1,:,ii)), squeeze(p2(2,:,ii)), squeeze(p2(3,:,ii)),'k'); 
    end
  else
    % plain surface
    for ii = 1:numel(knt1)
      plot (squeeze(p1(1,ii,:)), squeeze (p1(2,ii,:)),'k'); 
    end
    for ii = 1:numel(knt2)
      plot (p2(1,:,ii),p2(2,:,ii),'k');
    end
  end



 elseif (size (nurbs.knots,2) == 3) % plot a NURBS volume
   if (nargin < 2)
     nsub = [25 25 25];
   elseif (numel(nsub) == 1)
     nsub = [nsub nsub nsub];
   end
   % Plot the boundaries
   bnd = nrbextract (nurbs);
   nrbkntplot (bnd(1), nsub(2:3));
   hold on
   for iface = 2:6
     inds = setdiff(1:3, ceil(iface/2));
     nrbkntplot (bnd(iface), nsub(inds));
   end
 end
else % plot a NURBS curve
  if (nargin < 2)
    nsub = 1000;
  end
  nrbplot (nurbs, nsub);
  hold on

  % And plot the knots
   order = nurbs.order;
   p = nrbeval (nurbs, unique (nurbs.knots(order:end-order+1)));

   if (any (nurbs.coefs(3,:))) % plot a 3D curve
     plot3 (p(1,:), p(2,:), p(3,:), 'rx'); 
   else                     % plot a 2D curve
     plot (p(1,:), p(2,:), 'rx'); 
   end

end

if (~hold_flag)
  hold off
end

end
%!demo
%! crv = nrbtestcrv;
%! nrbkntplot(crv)
%! title('Test curve')
%! hold off

%!demo
%! sphere = nrbrevolve(nrbcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbkntplot(sphere);
%! title('Ball and torus - surface construction by revolution');
%! hold on;
%! torus = nrbrevolve(nrbcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbkntplot(torus);
%! hold off

%!demo
%! knots = {[0 0 0 1/2 1 1 1] [0 0 0 1 1 1]...
%!          [0 0 0 1/6 2/6 1/2 1/2 4/6 5/6 1 1 1]};
%!
%! coefs = [-1.0000   -0.9734   -0.7071    1.4290    1.0000    3.4172
%!          0    2.4172         0    0.0148   -2.0000   -1.9734
%!          0    2.0000    4.9623    9.4508    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -0.8536         0   -0.6036    1.9571    1.2071    3.5000
%!     0.3536    2.5000    0.2500    0.5429   -1.7071   -1.0000
%!          0    2.0000    4.4900    8.5444    3.4142    2.0000
%!     0.8536    1.0000    0.6036    1.0000    0.8536    1.0000
%!    -0.3536   -4.0000   -0.2500   -1.2929    1.7071    1.0000
%!     0.8536         0    0.6036   -2.7071   -1.2071   -5.0000
%!          0    2.0000    4.4900   10.0711    3.4142    2.0000
%!     0.8536    1.0000    0.6036    1.0000    0.8536    1.0000
%!          0   -4.0000         0    0.7071    2.0000    5.0000
%!     1.0000    4.0000    0.7071   -0.7071   -1.0000   -5.0000
%!          0    2.0000    4.9623   14.4142    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -2.5000   -4.0000   -1.7678    0.7071    1.0000    5.0000
%!          0    4.0000         0   -0.7071   -3.5000   -5.0000
%!          0    2.0000    6.0418   14.4142    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -2.4379         0   -1.7238    2.7071    1.9527    5.0000
%!     0.9527    4.0000    0.6737    1.2929   -3.4379   -1.0000
%!          0    2.0000    6.6827   10.0711    4.0000    2.0000
%!     1.0000    1.0000    0.7071    1.0000    1.0000    1.0000
%!    -0.9734   -1.0000   -0.6883    0.7071    3.4172    1.0000
%!     2.4172         0    1.7092   -1.4142   -1.9734   -2.0000
%!          0    4.0000    6.6827    4.9623    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!          0   -0.8536         0    0.8536    3.5000    1.2071
%!     2.5000    0.3536    1.7678   -1.2071   -1.0000   -1.7071
%!          0    3.4142    6.0418    4.4900    4.0000         0
%!     1.0000    0.8536    0.7071    0.6036    1.0000    0.8536
%!    -4.0000   -0.3536   -2.8284    1.2071    1.0000    1.7071
%!          0    0.8536         0   -0.8536   -5.0000   -1.2071
%!          0    3.4142    7.1213    4.4900    4.0000         0
%!     1.0000    0.8536    0.7071    0.6036    1.0000    0.8536
%!    -4.0000         0   -2.8284    1.4142    5.0000    2.0000
%!     4.0000    1.0000    2.8284   -0.7071   -5.0000   -1.0000
%!          0    4.0000   10.1924    4.9623    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!    -4.0000   -2.5000   -2.8284    0.7071    5.0000    1.0000
%!     4.0000         0    2.8284   -2.4749   -5.0000   -3.5000
%!          0    4.0000   10.1924    6.0418    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!          0   -2.4379         0    1.3808    5.0000    1.9527
%!     4.0000    0.9527    2.8284   -2.4309   -1.0000   -3.4379
%!          0    4.0000    7.1213    6.6827    4.0000         0
%!     1.0000    1.0000    0.7071    0.7071    1.0000    1.0000
%!    -1.0000   -0.9734    0.2071    2.4163    1.0000    3.4172
%!          0    2.4172   -1.2071   -1.3954   -2.0000   -1.9734
%!     2.0000    4.0000    7.0178    6.6827    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -0.8536         0    0.3536    2.4749    1.2071    3.5000
%!     0.3536    2.5000   -0.8536   -0.7071   -1.7071   -1.0000
%!     1.7071    4.0000    6.3498    6.0418    1.7071         0
%!     0.8536    1.0000    0.8536    0.7071    0.8536    1.0000
%!    -0.3536   -4.0000    0.8536    0.7071    1.7071    1.0000
%!     0.8536         0   -0.3536   -3.5355   -1.2071   -5.0000
%!     1.7071    4.0000    6.3498    7.1213    1.7071         0
%!     0.8536    1.0000    0.8536    0.7071    0.8536    1.0000
%!          0   -4.0000    1.2071    3.5355    2.0000    5.0000
%!     1.0000    4.0000   -0.2071   -3.5355   -1.0000   -5.0000
%!     2.0000    4.0000    7.0178   10.1924    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -2.5000   -4.0000   -0.5429    3.5355    1.0000    5.0000
%!          0    4.0000   -1.9571   -3.5355   -3.5000   -5.0000
%!     2.0000    4.0000    8.5444   10.1924    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000
%!    -2.4379         0   -0.0355    3.5355    1.9527    5.0000
%!     0.9527    4.0000   -1.4497   -0.7071   -3.4379   -1.0000
%!     2.0000    4.0000    9.4508    7.1213    2.0000         0
%!     1.0000    1.0000    1.0000    0.7071    1.0000    1.0000];
%! coefs = reshape (coefs, 4, 4, 3, 9);
%! horseshoe = nrbmak (coefs, knots);
%! nrbkntplot (horseshoe);
