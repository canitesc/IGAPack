function nrb = nrbreverse(nrb, idir)
%
% NRBREVERSE: Reverse the evaluation directions of a NURBS geometry.
% 
% Calling Sequence:
% 
%   rnrb = nrbreverse(nrb);
%   rnrb = nrbreverse(nrb, idir);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%   idir        : vector of directions to reverse.
%
% OUTPUT:
% 
%   rnrb	: Reversed NURBS.
% 
% Description:
% 
%   Utility function to reverse the evaluation direction of a NURBS
%   curve or surface.
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2013 Rafael Vazquez
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

if (nargin > 2)
  error('Incorrect number of input arguments');
end

if (iscell(nrb.knots))
  % reverse a NURBS surface or volume
  ndim = numel (nrb.knots);
  if (nargin == 1 || isempty (idir))
    idir = 1:ndim;
  end
  for ii = idir
    nrb.knots{ii} = sort (nrb.knots{ii}(end) - nrb.knots{ii});
    nrb.coefs = flip (nrb.coefs, ii+1);
  end

else
  % reverse a NURBS curve
  nrb.knots = sort (nrb.knots(end) - nrb.knots);
  nrb.coefs = fliplr (nrb.coefs);
end

end

%!demo
%! pnts = [0.5 1.5 3.0 7.5 8.5;
%!         3.0 5.5 1.5 4.0 4.5;
%!         0.0 0.0 0.0 0.0 0.0];
%! crv1 = nrbmak(pnts,[0 0 0 1/2 3/4 1 1 1]);
%! crv2 = nrbreverse(crv1);
%! fprintf('Knots of the original curve\n')
%! disp(crv1.knots)
%! fprintf('Knots of the reversed curve\n')
%! disp(crv2.knots)
%! fprintf('Control points of the original curve\n')
%! disp(crv1.coefs(1:2,:))
%! fprintf('Control points of the reversed curve\n')
%! disp(crv2.coefs(1:2,:))
%! nrbplot(crv1,100)
%! hold on
%! nrbplot(crv2,100)
%! title('The curve and its reverse are the same')
%! hold off

%!test
%! srf  = nrbrevolve(nrbline([1 0],[2 0]), [0 0 0], [0 0 1], pi/2);
%! srf  = nrbkntins (srf, {0.3, 0.6});
%! srf2 = nrbreverse (srf);
%! assert (srf.knots, cellfun(@(x) sort(1-x), srf2.knots, 'UniformOutput', false), 1e-15)
%! assert (srf.coefs, srf2.coefs(:,end:-1:1,end:-1:1))

%!test
%! srf  = nrbrevolve(nrbline([1 0],[2 0]), [0 0 0], [0 0 1], pi/2);
%! srf  = nrbkntins (srf, {0.3, 0.6});
%! srf2 = nrbreverse (srf, 1);
%! knt{1} = sort(1-srf2.knots{1}); knt{2} = srf2.knots{2};
%! assert (srf.knots, knt, 1e-15)
%! assert (srf.coefs, srf2.coefs(:,end:-1:1,:))
%! srf2 = nrbreverse (srf, 2);
%! knt{1} = srf2.knots{1}; knt{2} = sort(1-srf2.knots{2});
%! assert (srf.knots, knt, 1e-15)
%! assert (srf.coefs, srf2.coefs(:,:,end:-1:1))
