function varargout = nrbderiv (nurbs)
% 
% NRBDERIV: Construct up to the fourth derivative representation of a
%           NURBS curve, surface or volume.
% 
% Calling Sequence:
% 
%   ders = nrbderiv (nrb);
%   [ders, ders2] = nrbderiv (nrb);
%   [ders, ders2, ders3] = nrbderiv (nrb);
%   [ders, ders2, ders, ders4] = nrbderiv (nrb);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%
% OUTPUT:
% 
%   ders:  A data structure that represents the first
% 		    derivatives of a NURBS curve, surface or volume.
%   ders2: A data structure that represents the second
% 		    derivatives of a NURBS curve, surface or volume.
%   ders3: A data structure that represents the third
% 		    derivatives of a NURBS curve, surface or volume. (only surface
% 		    for now)
%   ders4: A data structure that represents the fourth
% 		    derivatives of a NURBS curve, surface or volume. (only surface
% 		    for now)
%
% 
% Description:
% 
%   The derivatives of a B-Spline are themselves a B-Spline of lower degree,
%   giving an efficient means of evaluating multiple derivatives. However,
%   although the same approach can be applied to NURBS, the situation for
%   NURBS is more complex. We have followed in this function the same idea
%   that was already used for the first derivative in the function nrbderiv.
%   The second derivative data structure can be evaluated later with the
%   function nrbdeval.
% 
% See also:
% 
%       nrbdeval
%
% Copyright (C) 2000 Mark Spink
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
% Copyright (C) 2018 Luca Coradello
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

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

% % % We raise the degree to avoid errors in the computation of the second
% % % derivative
% if (iscell (nurbs.knots))
%   ndim = size(nurbs.knots, 2);
% else
%   ndim = 1;
% end
% 
% if (nargout == 2)
%   degelev  = max (2*ones(1, ndim) - (nurbs.order-1), 0);
%   nurbs    = nrbdegelev (nurbs, degelev);
% end
% 
% degree = nurbs.order - 1;
% 
% if (ndim == 3)
% % NURBS structure represents a volume
%   num1 = nurbs.number(1);
%   num2 = nurbs.number(2);
%   num3 = nurbs.number(3);
% 
% % taking derivatives along the u direction
%   dknots = nurbs.knots;
%   dcoefs = permute (nurbs.coefs,[1 3 4 2]);
%   dcoefs = reshape (dcoefs,4*num2*num3,num1);
%   [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
%   dcoefs = permute (reshape (dcoefs,[4 num2 num3 size(dcoefs,2)]),[1 4 2 3]);
%   dnurbs{1} = nrbmak (dcoefs, dknots);
% 
%   if (nargout == 2)
% % taking second derivative along the u direction (duu)
%     dknots2 = dknots;
%     dcoefs2 = permute (dcoefs, [1 3 4 2]);
%     dcoefs2 = reshape (dcoefs2, 4*num2*num3, []);
%     [dcoefs2, dknots2{1}] = bspderiv (degree(1)-1, dcoefs2, dknots{1});
%     dcoefs2 = permute (reshape (dcoefs2, 4, num2, num3, []), [1 4 2 3]);
%     dnurbs2{1,1} = nrbmak (dcoefs2, dknots2); 
% 
% % taking second derivative along the v direction (duv and dvu)
%     dknots2 = dknots;
%     dcoefs2 = permute (dcoefs,[1 2 4 3]);
%     dcoefs2 = reshape (dcoefs2, 4*(num1-1)*num3, num2);
%     [dcoefs2, dknots2{2}] = bspderiv (degree(2), dcoefs2, dknots{2});
%     dcoefs2 = permute (reshape (dcoefs2, 4, num1-1, num3, []), [1 2 4 3]);
%     dnurbs2{1,2} = nrbmak (dcoefs2, dknots2);
%     dnurbs2{2,1} = dnurbs2{1,2};
% 
% % taking second derivative along the w direction (duw and dwu)
%     dknots2 = dknots;
%     dcoefs2 = reshape (dcoefs, 4*(num1-1)*num2, num3);
%     [dcoefs2, dknots2{3}] = bspderiv (degree(3), dcoefs2, dknots{3});
%     dcoefs2 = reshape (dcoefs2, 4, num1-1, num2, []);
%     dnurbs2{1,3} = nrbmak (dcoefs2, dknots2);
%     dnurbs2{3,1} = dnurbs2{1,3};
%   end
% 
% % taking derivatives along the v direction
%   dknots = nurbs.knots;
%   dcoefs = permute (nurbs.coefs,[1 2 4 3]);
%   dcoefs = reshape (dcoefs,4*num1*num3,num2);
%   [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
%   dcoefs = permute (reshape (dcoefs,[4 num1 num3 size(dcoefs,2)]),[1 2 4 3]);
%   dnurbs{2} = nrbmak (dcoefs, dknots);
% 
%   if (nargout == 2)
% % taking second derivative along the v direction (dvv)
%     dknots2 = dknots;
%     dcoefs2 = permute (dcoefs,[1 2 4 3]);
%     dcoefs2 = reshape (dcoefs2, 4*num1*num3, num2-1);
%     [dcoefs2, dknots2{2}] = bspderiv (degree(2)-1, dcoefs2, dknots{2});
%     dcoefs2 = permute (reshape (dcoefs2, 4, num1, num3, []), [1 2 4 3]);
%     dnurbs2{2,2} = nrbmak (dcoefs2, dknots2);
% 
% % taking second derivative along the w direction (dvw and dwv)
%     dknots2 = dknots;
%     dcoefs2 = reshape (dcoefs, 4*num1*(num2-1), num3);
%     [dcoefs2, dknots2{3}] = bspderiv (degree(3), dcoefs2, dknots{3});
%     dcoefs2 = reshape (dcoefs2, 4, num1, num2-1, []);
%     dnurbs2{2,3} = nrbmak (dcoefs2, dknots2);
%     dnurbs2{3,2} = dnurbs2{2,3};
%   end
% 
% % taking derivatives along the w direction
%   dknots = nurbs.knots;
%   dcoefs = reshape (nurbs.coefs,4*num1*num2,num3);
%   [dcoefs,dknots{3}] = bspderiv (degree(3),dcoefs,nurbs.knots{3});
%   dcoefs = reshape (dcoefs,[4 num1 num2 size(dcoefs,2)]);
%   dnurbs{3} = nrbmak (dcoefs, dknots);
% 
%   if (nargout == 2)
% % taking second derivative along the w direction (dww)
%     dknots2 = dknots;
%     dcoefs2 = reshape (dcoefs, 4*num1*num2, num3-1);
%     [dcoefs2, dknots2{3}] = bspderiv (degree(3)-1, dcoefs2, dknots{3});
%     dcoefs2 = reshape (dcoefs2, 4, num1, num2, []);
%     dnurbs2{3,3} = nrbmak (dcoefs2, dknots2);
%   end
% 
% elseif (ndim == 2)
% % NURBS structure represents a surface
%   num1 = nurbs.number(1);
%   num2 = nurbs.number(2);
% 
% % taking first derivative along the u direction
%   dknots = nurbs.knots;
%   dcoefs = permute (nurbs.coefs,[1 3 2]);
%   dcoefs = reshape (dcoefs,4*num2,num1);
%   [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
%   dcoefs = permute (reshape (dcoefs,[4 num2 size(dcoefs,2)]),[1 3 2]);
%   dnurbs{1} = nrbmak (dcoefs, dknots);
% 
%   if (nargout >= 2)
% % taking second derivative along the u direction (duu)
%     dknots2 = dknots;
%     dcoefs2 = permute (dcoefs, [1 3 2]);
%     dcoefs2 = reshape (dcoefs2, 4*num2, []);
%     [dcoefs2, dknots2{1}] = bspderiv (degree(1)-1, dcoefs2, dknots{1});
%     dcoefs2 = permute (reshape (dcoefs2, 4, num2, []), [1 3 2]);
%     dnurbs2{1,1} = nrbmak (dcoefs2, dknots2); 
% 
% % taking second derivative along the v direction (duv and dvu)
%     dknots2 = dknots;
%     dcoefs2 = reshape (dcoefs, 4*(num1-1), num2);
%     [dcoefs2, dknots2{2}] = bspderiv (degree(2), dcoefs2, dknots{2});
%     dcoefs2 = reshape (dcoefs2, 4, num1-1, []);
%     dnurbs2{1,2} = nrbmak (dcoefs2, dknots2);
%     dnurbs2{2,1} = dnurbs2{1,2};
%   end
% 
% % taking first derivative along the v direction
%   dknots = nurbs.knots;
%   dcoefs = reshape (nurbs.coefs,4*num1,num2);
%   [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
%   dcoefs = reshape (dcoefs,[4 num1 size(dcoefs,2)]);
%   dnurbs{2} = nrbmak (dcoefs, dknots);
% 
%   if (nargout >= 2)
% % taking second derivative along the v direction (dvv)
%     dknots2 = dknots;
%     dcoefs2 = reshape (dcoefs, 4*num1, num2-1);
%     [dcoefs2, dknots2{2}] = bspderiv (degree(2)-1, dcoefs2, dknots{2});
%     dcoefs2 = reshape (dcoefs2, 4, num1, []);
%     dnurbs2{2,2} = nrbmak (dcoefs2, dknots2);
%   end
% 
% else
%   % NURBS structure represents a curve
%   [dcoefs,dknots] = bspderiv (degree, nurbs.coefs, nurbs.knots);
%   dnurbs = nrbmak (dcoefs, dknots);
%   if (nargout == 2)
%     [dcoefs2,dknots2] = bspderiv (degree-1, dcoefs, dknots);
%     dnurbs2 = nrbmak (dcoefs2, dknots2);
%   end
% end


if (iscell (nurbs.knots))
  ndim = size(nurbs.knots, 2);
else
  ndim = 1;
end
% We raise the degree to avoid errors in the computation of the higher
% order derivatives
if (nargout >= 2)
  degelev  = max (nargout*ones(1, ndim) - (nurbs.order-1), 0);
  nurbs    = nrbdegelev (nurbs, degelev);
end

if (ndim == 1) % in case of a curve create a cell array to use a dimension-indipendent algorithm
  tmp = nurbs.knots;
  nurbs.knots = {};
  nurbs.knots{1} = tmp;
end

for idim = 1:ndim
    num = nurbs.number;
    degree = nurbs.order - 1;
    dknots = nurbs.knots;
    coord = idim + 1;
    vec_permute = [setdiff(1:ndim+1,coord),coord];
    dcoefs = permute (nurbs.coefs, vec_permute);
    vec_reshape = setdiff(1:ndim,idim);
    dcoefs = reshape (dcoefs,4*prod(num(vec_reshape)), num(idim));
    [dcoefs, dknots{idim}] = bspderiv (degree(idim), dcoefs, nurbs.knots{idim});
    dcoefs = reshape (dcoefs, [4 num(vec_reshape) num(idim)-1]);
    [~,~,ib] = intersect(1:ndim+1, vec_permute);
    dcoefs = permute(dcoefs, ib');
    dnurbs{idim} = nrbmak (dcoefs, dknots);

    if (nargout >= 2) % second derivatives
        degree(idim) = degree(idim) - 1;
        num(idim) = num(idim) - 1;
        for jdim = 1:ndim
            dknots2 = dknots;
            coord = jdim + 1;
            vec_permute = [setdiff(1:ndim+1,coord),coord];
            dcoefs2 = permute (dcoefs, vec_permute);
            vec_reshape = setdiff(1:ndim,jdim);
            dcoefs2 = reshape (dcoefs2, 4*prod(num(vec_reshape)), num(jdim));
            [dcoefs2, dknots2{jdim}] = bspderiv (degree(jdim), dcoefs2, dknots{jdim});
            dcoefs2 = reshape (dcoefs2, [4 num(vec_reshape) num(jdim)-1]);
            [~,~,ib] = intersect(1:ndim+1, vec_permute);
            dcoefs2 = permute(dcoefs2, ib');
            dnurbs2{idim,jdim} = nrbmak (dcoefs2, dknots2); 
            
            if(nargout >= 3) %third derivatives   
                degree(jdim) = degree(jdim) - 1;
                num(jdim) = num(jdim) - 1;
                for kdim = 1:ndim
                    dknots3 = dknots2;
                    coord = kdim + 1;
                    vec_permute = [setdiff(1:ndim+1,coord),coord];
                    dcoefs3 = permute (dcoefs2, vec_permute);
                    vec_reshape = setdiff(1:ndim,kdim);
                    dcoefs3 = reshape (dcoefs3, 4*prod(num(vec_reshape)), num(kdim));
                    [dcoefs3, dknots3{kdim}] = bspderiv (degree(kdim), dcoefs3, dknots2{kdim});
                    dcoefs3 = reshape (dcoefs3, [4 num(vec_reshape) num(kdim)-1]);
                    [~,~,ib] = intersect(1:ndim+1, vec_permute);
                    dcoefs3 = permute(dcoefs3, ib');
                    dnurbs3{idim,jdim,kdim} = nrbmak (dcoefs3, dknots3);
                    
                    if(nargout == 4) %fourth derivatives    
                        degree(kdim) = degree(kdim) - 1;
                        num(kdim) = num(kdim) - 1;
                        for ldim = 1:ndim
                            dknots4 = dknots3;
                            coord = ldim + 1;
                            vec_permute = [setdiff(1:ndim+1,coord),coord];
                            dcoefs4 = permute (dcoefs3, vec_permute);
                            vec_reshape = setdiff(1:ndim,ldim);
                            dcoefs4 = reshape (dcoefs4, 4*prod(num(vec_reshape)), num(ldim));
                            [dcoefs4, dknots4{ldim}] = bspderiv (degree(ldim), dcoefs4, dknots3{ldim});
                            dcoefs4 = reshape (dcoefs4, [4 num(vec_reshape) num(ldim)-1]);
                            [~,~,ib] = intersect(1:ndim+1, vec_permute);
                            dcoefs4 = permute(dcoefs4, ib');
                            dnurbs4{idim,jdim,kdim,ldim} = nrbmak (dcoefs4, dknots4);
%                             degree(ldim) = degree(ldim) - 1;
%                             num(ldim) = num(ldim) - 1;
                        end 
                        degree(kdim) = degree(kdim) + 1;
                        num(kdim) = num(kdim) + 1;

                    end
                end
                degree(jdim) = degree(jdim) + 1;
                num(jdim) = num(jdim) + 1;
            end
        end
    end
end


if (ndim == 1) % in the case of a curve transform everything back otherwise it will throw errors somewhere else in the code !
    tmp = nurbs.knots{1};
    nurbs.knots = tmp;
    tmp = dnurbs{1};
    dnurbs = tmp;
    if (nargout >= 2)
        tmp = dnurbs2{1,1};
        dnurbs2 = tmp;
    end
    if (nargout >= 3)
        tmp = dnurbs3{1,1,1};
        dnurbs3 = tmp;
    end
    if (nargout == 4)
        tmp = dnurbs4{1,1,1,1};
        dnurbs4 = tmp;
    end

end


varargout{1} = dnurbs;
if (nargout >= 2)
  varargout{2} = dnurbs2;
  
  if (iscell (dnurbs2))
    dnurbs2 = [dnurbs2{:}];
  end
  if (any (arrayfun(@(x) any(isnan(x.coefs(:)) | isinf(x.coefs(:))), dnurbs2)))
    warning ('nrbderiv:SecondDerivative', ...
        ['The structure for the second derivative contains Inf and/or NaN coefficients, ' ...
        'probably due to low continuity at repeated knots. This should not affect the ' ...
        'computation of the second derivatives, except at those knots.'])
  end
end
if(nargout >= 3)
  varargout{3} = dnurbs3;
end
if(nargout == 4)
  varargout{4} = dnurbs4;
end
end

%!demo
%! crv = nrbtestcrv;
%! nrbplot(crv,48);
%! title('First derivatives along a test curve.');
%! 
%! tt = linspace(0.0,1.0,9);
%! 
%! dcrv = nrbderiv(crv);
%! 
%! [p1, dp] = nrbdeval(crv,dcrv,tt);
%! 
%! p2 = vecnormalize(dp);
%! 
%! hold on;
%! plot(p1(1,:),p1(2,:),'ro');
%! h = quiver(p1(1,:),p1(2,:),p2(1,:),p2(2,:),0);
%! set(h,'Color','black');
%! hold off;

%!demo
%! srf = nrbtestsrf;
%! p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)});
%! h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
%! set(h,'FaceColor','blue','EdgeColor','blue');
%! title('First derivatives over a test surface.');
%!
%! npts = 5;
%! tt = linspace(0.0,1.0,npts);
%! dsrf = nrbderiv(srf);
%! 
%! [p1, dp] = nrbdeval(srf, dsrf, {tt, tt});
%! 
%! up2 = vecnormalize(dp{1});
%! vp2 = vecnormalize(dp{2});
%! 
%! hold on;
%! plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
%! h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:));
%! h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:));
%! set(h1,'Color','black');
%! set(h2,'Color','black');
%! 
%! hold off;

%!test
%! knots = [0 0 0 0.5 1 1 1];
%! coefs(1,:) = [0 2 4 2];
%! coefs(2,:) = [0 2 2 0];
%! coefs(3,:) = [0 4 2 0];
%! coefs(4,:) = [1 2 2 1];
%! nrb = nrbmak (coefs, knots);
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! x = linspace (0, 1, 10);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, x);
%! w = -4*x.^2 + 4*x + 1;
%! F = zeros (3,numel(x)); DF = zeros (3, numel(x)); D2F = zeros (3, numel(x));
%! F(1,:) = (-4*x.*(x-2)./w) .* (x<0.5) + ((4*x - 5)./w + 3) .* (x>0.5);
%! F(2,:) = (2-2./w);
%! F(3,:) = (-4*x.*(5*x-4)./w) .* (x<0.5) + (-4*(x.^2 - 1)./w) .* (x>0.5);
%! DF(1,:) = (8*(2*x.^2-x+1)./w.^2) .* (x<0.5) + (8*(2*x-3).*(x-1)./w.^2) .* (x>0.5);
%! DF(2,:) = -8*(2*x-1)./w.^2;
%! DF(3,:) = -(8*(2*x.^2+5*x-2)./w.^2) .* (x<0.5) - (8*(2*x.^2-3*x+2)./w.^2) .* (x>0.5);
%! D2F(1,:) = 8*(16*x.^3-12*x.^2+24*x-9)./w.^3 .* (x<0.5) + 8*(16*x.^3-60*x.^2+72*x-29)./w.^3 .* (x>0.5);
%! D2F(2,:) = -16*(12*x.^2-12*x+5)./w.^3;
%! D2F(3,:) = -8*(16*x.^3+60*x.^2-48*x+21)./w.^3 .* (x<0.5) -8*(16*x.^3-36*x.^2+48*x-19)./w.^3 .* (x>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (DF, jac, 1e3*eps)
%! assert (D2F, hess, 1e3*eps)

%!test
%! knots = {[0 0 0 1 1 1], [0 0 0 0.5 1 1 1]};
%! coefs = ones (4,3,4);
%! coefs(1,:,:) = reshape ([0 0 0 0; 1 1 1 1; 2 2 4 2], 1, 3, 4);
%! coefs(2,:,:) = reshape ([0 1 2 3; 0 1 2 3; 0 1 4 3], 1, 3, 4);
%! coefs(3,:,:) = reshape ([0 1 0 0; 0 0 0 0; 0 0 0 0], 1, 3, 4);
%! coefs(4,:,:) = reshape ([1 1 1 1; 1 1 1 1; 1 1 2 1], 1, 3, 4);
%! nrb = nrbmak (coefs, knots);
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 4); Y = linspace (0, 1, 4);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y});
%! [y, x] = meshgrid (X, Y);
%! w = (2*x.^2.*y.^2 + 1) .* (y < 0.5) + (-6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 1) .* (y > 0.5);
%! F = zeros ([3,size(x)]);
%! F(1,:,:) = ((2*x - 2) ./w + 2) .* (y<0.5) + (2 + (2*x-2)./w) .* (y > 0.5);
%! F(2,:,:) = (2 - (2*(y-1).^2)./w).*(y<0.5) + ...
%!     ((-12*x.^2.*y.^2 + 16*x.^2.*y - 4*x.^2 + 2*y.^2 + 1)./w).*(y>0.5);
%! F(3,:,:) = (-2*y.*(3*y - 2).*(x - 1).^2./w) .* (y<0.5) + ...
%!     (2*(x - 1).^2.*(y - 1).^2./w) .* (y>0.5); 
%! dFdu = zeros ([3,size(x)]);
%! dFdu(1,:,:) = (((8*x - 4*x.^2).*y.^2 + 2)./w.^2).*(y<0.5) + ...
%!     (((12*y.^2 - 16*y + 4).*x.^2 + (-24*y.^2 + 32*y - 8).*x + 2)./w.^2).*(y>0.5);
%! dFdu(2,:,:) = (8*x.*y.^2.*(y - 1).^2./w.^2).*(y<0.5) + ...
%!     ((4*x.*(3*y - 1).*(2*y.^2 - 1).*(y - 1))./w.^2).*(y>0.5);
%! dFdu(3,:,:) = (-4*y.*(2.*x.*y.^2 + 1).*(3*y - 2).*(x - 1)./w.^2).*(y<0.5) + ...
%!     ((-4*(x - 1).*(y - 1).^2.*(6*x.*y.^2 - 8*x.*y + 2*x - 1))./w.^2).*(y>0.5);
%! dFdv = zeros ([3,size(x)]);
%! dFdv(1,:,:) = (-8*x.^2.*y.*(x - 1)./w.^2).*(y<0.5) + ...
%!     (8*x.^2.*(3*y - 2).*(x - 1)./w.^2).*(y>0.5);
%! dFdv(2,:,:) = (-4*(2*y.*x.^2 + 1).*(y - 1)./w.^2).*(y<0.5) + ...
%!     (((16*y.^2 - 20*y + 8).*x.^2 + 4*y)./w.^2).*(y>0.5);
%! dFdv(3,:,:) = (-4*(x - 1).^2.*(2*x.^2.*y.^2 + 3*y - 1)./w.^2).*(y<0.5) + ...
%!     (4*(x - 1).^2.*(y - 1).*(2*x.^2 - 2*x.^2.*y + 1)./w.^2).*(y>0.5);
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:) = (-((48*x.^2 - 16*x.^3).*y.^4 + (24*x - 8).*y.^2)./w.^3).*(y<0.5) + ...
%!     (((32*(3*y - 1).*(x - 1).*(y - 1))-(8*(3*y - 1).*(x - 3).*(y - 1).*w))./w.^3).*(y>0.5);
%! d2Fduu(2,:,:) = (-(8*y.^2.*(6*x.^2.*y.^2 - 1).*(y - 1).^2)./w.^3).*(y<0.5) + ...
%!     ((4*(3*y - 1).*(2*y.^2 - 1).*(y - 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 + 1))./w.^3).*(y>0.5);
%! d2Fduu(3,:,:) =  ((4*y.*(3*y - 2).*(8*x.^3.*y.^4 - 12*x.^2.*y.^4 + 6*x.^2.*y.^2 - 12*x.*y.^2 + 2*y.^2 - 1))./w.^3).*(y<0.5) + ...
%!  ((4*(y - 1).^2.*(6*y.^2 - 8*y + 3) - 4*x.^3.*(y - 1).^2.*(72*y.^4 - 192*y.^3 + 176*y.^2 - 64*y + 8) + 4*x.^2.*(y - 1).^2.*(108*y.^4 - 288*y.^3 + 282*y.^2 - 120*y + 18) - 4*x.*(y - 1).^2.*(36*y.^2 - 48*y + 12))./w.^3) .* (y>0.5);
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:) = (8*x.^2.*(6*x.^2.*y.^2 - 1).*(x - 1)./w.^3) .* (y<0.5) + ...
%!      (8*x.^2.*(x - 1).*(54*x.^2.*y.^2 - 72*x.^2.*y + 26*x.^2 + 3)./w.^3) .* (y>0.5);
%! d2Fdvv(2,:,:) =  (-((48*y.^2 - 32*y.^3).*x.^4 + (- 24*y.^2 + 48*y - 8).*x.^2 + 4)./w.^3) .*(y<0.5) + ...
%!       (((192*y.^3 - 360*y.^2 + 288*y - 88).*x.^4 + (72*y.^2 - 28).*x.^2 + 4)./w.^3) .* (y>0.5);
%! d2Fdvv(3,:,:) =  (4*(x - 1).^2.*(8*x.^4.*y.^3 + 18*x.^2.*y.^2 - 12*x.^2.*y - 3))./w.^3 .* (y<0.5) + ...
%!     ((4*(x - 1).^2.*(24*x.^4 + 18*x.^2 + 1) + 4*y.^2.*(72*x.^4 + 18*x.^2).*(x - 1).^2 - 96*x.^4.*y.^3.*(x - 1).^2 - 4*y.*(72*x.^4 + 36*x.^2).*(x - 1).^2)./w.^3) .* (y>0.5);
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:) = (-(y.^3.*(32*x.^3 - 16*x.^4) - y.*(16*x - 24*x.^2))./w.^3) .* (y<0.5) + ...
%!     (-(-8*(3*y - 2).*(6*y.^2 - 8*y + 2).*x.^4 + 8*(3*y - 2).*(12*y.^2 - 16*y + 4).*x.^3 + (48 - 72*y).*x.^2 + (48*y - 32).*x)./w.^3) .* (y>0.5);
%! d2Fduv(2,:,:) = (16*x.*y.*(y - 1).*(2*x.^2.*y.^2 + 2*y - 1)./w.^3) .* (y<0.5) + ...
%!     (-(8*x.*(4*y.^2 - 5*y + 2))./w.^2 + (16*x.*(3*y - 2).*(2*y.^2 - 1))./w.^3) .* (y>0.5);
%! d2Fduv(3,:,:) = (-(8*(x - 1).*(4*x.^3.*y.^4 - 6*x.^2.*y.^3 + 6*x.^2.*y.^2 + 12*x.*y.^3 - 6*x.*y.^2 + 3*y - 1))./w.^3) .* (y<0.5) + ...
%!     ((8*(x - 1).*(y - 1).*(12*x.^3.*y.^3 - 28*x.^3.*y.^2 + 20*x.^3.*y - 4*x.^3 + 6*x.^2.*y.^2 - 12*x.^2.*y + 6*x.^2 - 12*x.*y.^2 + 18*x.*y - 6*x + 1))./w.^3) .* (y>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (dFdu, jac{1}, 1e3*eps)
%! assert (dFdv, jac{2}, 1e3*eps)
%! assert (d2Fduu, hess{1,1}, 1e3*eps)
%! assert (d2Fduv, hess{1,2}, 1e3*eps)
%! assert (d2Fduv, hess{2,1}, 1e3*eps)
%! assert (d2Fdvv, hess{2,2}, 1e3*eps)

%!test
%! knots = {[0 0 0 1 1 1], [0 0 0 0.5 1 1 1]};
%! coefs = ones (4,3,4);
%! coefs(1,:,:) = reshape ([0 0 0 0; 1 1 1 1; 2 2 4 2], 1, 3, 4);
%! coefs(2,:,:) = reshape ([0 1 2 3; 0 1 2 3; 0 1 4 3], 1, 3, 4);
%! coefs(3,:,:) = reshape ([0 1 0 0; 0 0 0 0; 0 0 0 0], 1, 3, 4);
%! coefs(4,:,:) = reshape ([1 1 1 1; 1 1 1 1; 1 1 2 1], 1, 3, 4);
%! nrb = nrbmak (coefs, knots);
%! nrb = nrbdegelev (nrbextrude (nrb, [0.4 0.6 2]), [0 0 1]);
%! nrb.coefs(4,2,3,3) = 1.5;
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 4); Y = linspace (0, 1, 4); Z = linspace (0, 1, 4);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y Z});
%! [y, x, z] = meshgrid (X, Y, Z);
%! w = (-2*x.^2.*y.^2.*z.^2 + 2*x.^2.*y.^2 + 2*x.*y.^2.*z.^2 + 1) .* (y < 0.5) + ...
%!     (6*x.^2.*y.^2.*z.^2 - 6*x.^2.*y.^2 - 8*x.^2.*y.*z.^2 + 8*x.^2.*y + 2*x.^2.*z.^2 - 2*x.^2 - 6*x.*y.^2.*z.^2 + 8*x.*y.*z.^2 - 2*x.*z.^2 + 1) .* (y > 0.5);
%! F = zeros ([3,size(x)]);
%! F(1,:,:,:) = ((10*x + 20*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2))./(5*w)) .* (y<0.5) + ...
%!     (60*x.^2.*y.^2 - 10*x + z.*(12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2) - 80*x.^2.*y + 20*x.^2)./(-5*w) .* (y > 0.5);
%! F(2,:,:,:) = ((20*y + 20*x.^2.*y.^2 + z.*(6*x.^2.*y.^2 + 3) - 10*y.^2)./(5*w)).*(y<0.5) + ...
%!     ((60*x.^2.*y.^2 + z.*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3) - 80*x.^2.*y + 20*x.^2 - 10*y.^2 - 5)./(-5*w)).*(y>0.5);
%! F(3,:,:,:) = ((4*y - 6*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2) - 8*x.*y + 12*x.*y.^2 + 4*x.^2.*y - 6*y.^2)./w) .* (y<0.5) + ...
%!     ((2*z - 4*y - 4*x + 2*x.^2.*y.^2 + 8*x.*y - 4*x.*y.^2 - 4*x.^2.*y - 4*x.^2.*z + 2*x.^2 + 2*y.^2 + 16*x.^2.*y.*z - 12*x.^2.*y.^2.*z + 2)./w) .* (y>0.5); 
%! dFdu = zeros ([3,size(x)]);
%! dFdu(1,:,:,:) = ((x.*((8*y.^2.*z.^3)/5 + 8*y.^2) - (4*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 2)./w.^2).*(y<0.5) + ...
%!     ((z.^3.*(x.^2.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (16*y)/5 - x.*((24*y.^2)/5 - (32*y)/5 + 8/5) + (12*y.^2)/5 + 4/5) - x.*(24*y.^2 - 32*y + 8) + x.^2.*(12*y.^2 - 16*y + 4) + x.^2.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) + 2)./w.^2).*(y>0.5);
%! dFdu(2,:,:,:) = ((z.^2.*(8*x.^2.*y.^4 - y.^2.*(8*y - 4*y.^2) + (2*x.*y.^2.*(40*y - 20*y.^2))/5) + z.^3.*((12*x.^2.*y.^4)/5 + (12*x.*y.^2)/5 - (6*y.^2)/5) + (2*x.*y.^2.*(20*y.^2 - 40*y + 20))/5)./w.^2).*(y<0.5) + ...
%!     (((2*(3*y.^2 - 4*y + 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 6*x + 3).*z.^3)/5 + (2*(3*y.^2 - 4*y + 1).*(60*x.^2.*y.^2 - 80*x.^2.*y + 20*x.^2 - 20*x.*y.^2 - 10*x + 10*y.^2 + 5).*z.^2)/5 - (2*(10*x - 20*x.*y.^2).*(3*y.^2 - 4*y + 1))/5)./w.^2).*(y>0.5);
%! dFdu(3,:,:,:) = ((4*y.*(3*y - 2) + z.^3.*(8*x.^2.*y.^4 + 8*x.*y.^2 - 4*y.^2) - z.^2.*(4*y.*(2*y.^2 - 3*y.^3).*x.^2 - 4*y.*(4*y.^2 - 6*y.^3).*x + 4*y.*(2*y.^2 - 3*y.^3)) + 4*x.^2.*y.*(4*y.^2 - 6*y.^3) - 4*x.*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2)) ./w.^2).*(y<0.5) + ...
%!     ((z.^2.*(4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1).*x.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2).*x + 4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)) - 4*(y - 1).^2 + z.^3.*(4*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2).*x.^2 - 4*(6*y - 2).*(y - 1).*x + 4*(3*y - 1).*(y - 1)) + 4*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) - 4*x.^2.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2))./w.^2) .* (y > 0.5); 
%! dFdv = zeros ([3,size(x)]);
%! dFdv(1,:,:,:) = ((8*x.*y.*(x - 1).*(z.^3 + 5*x.*z.^2 - 5*x))/5./w.^2).*(y<0.5) + ...
%!     (-(8*x.*(3*y - 2).*(x - 1).*(z.^3 + 5*x.*z.^2 - 5*x))/5./w.^2).*(y>0.5);
%! dFdv(2,:,:,:) = (-((8*x.*z.^2 - x.^2.*(8*z.^2 - 8)).*y.^2 + ((12*x.*z.^3)/5 - x.^2.*((12*z.^3)/5 + 8) + 4).*y - 4)./w.^2).*(y<0.5) + ...
%!     ((4*y + z.^3.*(x.*((36*y)/5 - 24/5) - x.^2.*((36*y)/5 - 24/5)) + z.^2.*(x.*(16*y.^2 + 4*y - 8) - x.^2.*(16*y.^2 + 4*y - 8)) + x.^2.*(16*y.^2 - 20*y + 8))./w.^2).*(y>0.5);
%! dFdv(3,:,:,:) = ((4*(x - 1).^2 - y.*(4*(3*x - 3).*(x - 1) - 8*x.*z.^3.*(x - 1)) + y.^2.*(4*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x).*z.^2 + 4*(2*x.^2 - 2*x.^3).*(x - 1)))./w.^2).*(y<0.5) + ...
%!     ((y.^2.*(4*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x).*z.^2 + 4*(2*x.^2 - 2*x.^3).*(x - 1)) - 4*(x - 1).*(2*x.^3 - 2*x.^2 + x - 1) - y.*(24*x.*(x - 1).*z.^3 + 4*(x - 1).*(4*x.^3 - 8*x.^2 + 4*x).*z.^2 - 4*(x - 1).*(4*x.^3 - 4*x.^2 + x - 1)) + 16*x.*z.^3.*(x - 1) + 4*z.^2.*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x))./w.^2).*(y>0.5);
%! dFdw = zeros ([3,size(x)]);
%! dFdw(1,:,:,:) = ((4*x.^2.*y.^2 + 2)./(- 10*x.^2.*y.^2.*z.^2 + 10*x.^2.*y.^2 + 10*x.*y.^2.*z.^2 + 5) - ((20*x.*y.^2.*z - 20*x.^2.*y.^2.*z).*(10*x + 20*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2)))./(5*w).^2).*(y<0.5) + ...
%!     ((12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2)./(- 30*x.^2.*y.^2.*z.^2 + 30*x.^2.*y.^2 + 40*x.^2.*y.*z.^2 - 40*x.^2.*y - 10*x.^2.*z.^2 + 10*x.^2 + 30*x.*y.^2.*z.^2 - 40*x.*y.*z.^2 + 10*x.*z.^2 - 5) - ((60*x.^2.*y.^2 - 10*x + z.*(12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2) - 80*x.^2.*y + 20*x.^2).*(- 60*z.*x.^2.*y.^2 + 80*z.*x.^2.*y - 20*z.*x.^2 + 60*z.*x.*y.^2 - 80*z.*x.*y + 20*z.*x))./(5*w).^2).*(y>0.5);
%! dFdw(2,:,:,:) = ((6*x.^2.*y.^2 + 3)./(- 10*x.^2.*y.^2.*z.^2 + 10*x.^2.*y.^2 + 10*x.*y.^2.*z.^2 + 5) - ((20*x.*y.^2.*z - 20*x.^2.*y.^2.*z).*(20*y + 20*x.^2.*y.^2 + z.*(6*x.^2.*y.^2 + 3) - 10*y.^2))./(5*w).^2).*(y<0.5) + ...
%!     ((18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3)./(- 30*x.^2.*y.^2.*z.^2 + 30*x.^2.*y.^2 + 40*x.^2.*y.*z.^2 - 40*x.^2.*y - 10*x.^2.*z.^2 + 10*x.^2 + 30*x.*y.^2.*z.^2 - 40*x.*y.*z.^2 + 10*x.*z.^2 - 5) - ((- 60*z.*x.^2.*y.^2 + 80*z.*x.^2.*y - 20*z.*x.^2 + 60*z.*x.*y.^2 - 80*z.*x.*y + 20*z.*x).*(60*x.^2.*y.^2 + z.*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3) - 80*x.^2.*y + 20*x.^2 - 10*y.^2 - 5))./(5*w).^2).*(y>0.5);
%! dFdw(3,:,:,:) = ((4*x.^2.*y.^2 + 2)./(2*x.^2.*y.^2 - z.^2.*(2*x.^2.*y.^2 - 2*x.*y.^2) + 1) + (2*z.*(2*x.^2.*y.^2 - 2*x.*y.^2).*(4*y - 6*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2) - 8*x.*y + 12*x.*y.^2 + 4*x.^2.*y - 6*y.^2))./w.^2).*(y<0.5) + ...
%!     ((12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2)./(6*x.^2.*y.^2 + z.^2.*(- 6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 6*x.*y.^2 - 8*x.*y + 2*x) - 8*x.^2.*y + 2*x.^2 - 1) + (2*z.*(- 6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 6*x.*y.^2 - 8*x.*y + 2*x).*(2*z - 4*y - 4*x + 2*x.^2.*y.^2 + 8*x.*y - 4*x.*y.^2 - 4*x.^2.*y - 4*x.^2.*z + 2*x.^2 + 2*y.^2 + 16*x.^2.*y.*z - 12*x.^2.*y.^2.*z + 2))./w.^2).*(y>0.5);
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:,:) = (((8*y.^2.*z.^3)/5 + 2*x.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 8*y.^2)./w.^2 - (2*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2).*(x.*((8*y.^2.*z.^3)/5 + 8*y.^2) - (4*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 2))./w.^3).*(y<0.5) + ...
%!     ((32*y + 2*x.*(12*y.^2 - 16*y + 4) + z.^3.*((32*y)/5 + 2*x.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (24*y.^2)/5 - 8/5) - 24*y.^2 + 2*x.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) - 8)./w.^2 - (2*(z.^3.*(x.^2.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (16*y)/5 - x.*((24*y.^2)/5 - (32*y)/5 + 8/5) + (12*y.^2)/5 + 4/5) - x.*(24*y.^2 - 32*y + 8) + x.^2.*(12*y.^2 - 16*y + 4) + x.^2.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) + 2).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3).*(y>0.5);
%! d2Fduu(2,:,:,:) = ((z.^3.*((24*x.*y.^4)/5 + (12*y.^2)/5) + (2*y.^2.*(20*y.^2 - 40*y + 20))/5 + z.^2.*((2*y.^2.*(40*y - 20*y.^2))/5 + 16*x.*y.^4))./w.^2 - (2*(z.^2.*(8*x.^2.*y.^4 - y.^2.*(8*y - 4*y.^2) + (2*x.*y.^2.*(40*y - 20*y.^2))/5) + z.^3.*((12*x.^2.*y.^4)/5 + (12*x.*y.^2)/5 - (6*y.^2)/5) + (2*x.*y.^2.*(20*y.^2 - 40*y + 20))/5).*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2))./w.^3).*(y<0.5) + ...
%!     (((2*(3*y.^2 - 4*y + 1).*(36*x.*y.^2 - 48*x.*y + 12*x - 6).*z.^3)/5 - (2*(3*y.^2 - 4*y + 1).*(160*x.*y - 40*x - 120*x.*y.^2 + 20*y.^2 + 10).*z.^2)/5 + (2*(20*y.^2 - 10).*(3*y.^2 - 4*y + 1))/5)./w.^2 - (2*((2*(3*y.^2 - 4*y + 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 6*x + 3).*z.^3)/5 + (2*(3*y.^2 - 4*y + 1).*(60*x.^2.*y.^2 - 80*x.^2.*y + 20*x.^2 - 20*x.*y.^2 - 10*x + 10*y.^2 + 5).*z.^2)/5 - (2*(10*x - 20*x.*y.^2).*(3*y.^2 - 4*y + 1))/5).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3).*(y>0.5);
%! d2Fduu(3,:,:,:) =  (((16*x.*y.^4 + 8*y.^2).*z.^3 + (4*y.*(4*y.^2 - 6*y.^3) - 8*x.*y.*(2*y.^2 - 3*y.^3)).*z.^2 - 4*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2) + 8*x.*y.*(4*y.^2 - 6*y.^3))./w.^2 - (2*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2).*(4*y.*(3*y - 2) + z.^3.*(8*x.^2.*y.^4 + 8*x.*y.^2 - 4*y.^2) - z.^2.*(4*y.*(2*y.^2 - 3*y.^3).*x.^2 - 4*y.*(4*y.^2 - 6*y.^3).*x + 4*y.*(2*y.^2 - 3*y.^3)) + 4*x.^2.*y.*(4*y.^2 - 6*y.^3) - 4*x.*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2)))./w.^3).*(y<0.5) + ...
%!  (-((4*(6*y - 2).*(y - 1) - 8*x.*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2)).*z.^3 + (4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2) - 8*x.*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)).*z.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) + 8*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2))./w.^2 - (2*(z.^2.*(4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1).*x.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2).*x + 4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)) - 4*(y - 1).^2 + z.^3.*(4*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2).*x.^2 - 4*(6*y - 2).*(y - 1).*x + 4*(3*y - 1).*(y - 1)) + 4*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) - 4*x.^2.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2)).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:,:) = ((((8.*x.^2.*(6.*z.^3 - 6.*z.^5))/5 + (8.*x.^4.*(10.*z.^4 - 20.*z.^2 + 10))/5 - (8.*x.^3.*(- 4.*z.^5 + 10.*z.^4 + 4.*z.^3 - 30.*z.^2 + 20))/5 + (16.*x.*z.^5)/5).*y.^3 + ((8.*x.*(2.*z.^3 - 10.*z.^2 + 10))/5 + (8.*x.^2.*(15.*z.^2 - 15))/5 - (8.*z.^3)/5).*y)./w.^3) .* (y<0.5) + ...
%!     (-(x.^4.*((8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10).*z.^4)/5 - (8.*(3.*y - 2).*(60.*y.^2 - 80.*y + 20).*z.^2)/5 + (8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10))/5) - x.^3.*(- (8.*(3.*y - 2).*(12.*y.^2 - 16.*y + 4).*z.^5)/5 + (8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10).*z.^4)/5 + (8.*(3.*y - 2).*(12.*y.^2 - 16.*y + 4).*z.^3)/5 - (8.*(3.*y - 2).*(90.*y.^2 - 120.*y + 30).*z.^2)/5 + (8.*(3.*y - 2).*(60.*y.^2 - 80.*y + 20))/5) + z.^3.*((24.*y)/5 - 16/5) - x.^2.*((8.*(3.*y - 2).*(18.*y.^2 - 24.*y + 6).*z.^5)/5 - (8.*(3.*y - 2).*(18.*y.^2 - 24.*y + 6).*z.^3)/5 + (72.*y - 48).*z.^2 - 72.*y + 48) + x.*((8.*(3.*y - 2).*(6.*y.^2 - 8.*y + 2).*z.^5)/5 + (32/5 - (48.*y)/5).*z.^3 + (48.*y - 32).*z.^2 - 48.*y + 32))./(-w).^3) .* (y>0.5);
%! d2Fduv(2,:,:,:) = ((((4.*x.^2.*(60.*z.^2 - 60.*z.^4))/5 + (4.*x.^3.*(40.*z.^4 - 80.*z.^2 + 40))/5 + 16.*x.*z.^4).*y.^4 + ((4.*x.^2.*(18.*z.^3 - 18.*z.^5))/5 + (4.*x.^3.*(12.*z.^5 - 12.*z.^3 + 40.*z.^2 - 40))/5 + (4.*x.*(6.*z.^5 - 40.*z.^2 + 40))/5 + 16.*z.^2).*y.^3 + ((4.*x.*(60.*z.^2 - 60))/5 - 24.*z.^2).*y.^2 + ((4.*x.*(6.*z.^3 + 20))/5 - (12.*z.^3)/5).*y)./w.^3) .* (y<0.5) + ...
%!     ((z.^3.*(((432.*y.^3)/5 - (864.*y.^2)/5 + (528.*y)/5 - 96/5).*x.^3 + (- (648.*y.^3)/5 + (1296.*y.^2)/5 - (792.*y)/5 + 144/5).*x.^2 + ((72.*y)/5 - 48/5).*x - (36.*y)/5 + 24/5) - x.^3.*(192.*y.^4 - 496.*y.^3 + 480.*y.^2 - 208.*y + 32) + z.^4.*((- 192.*y.^4 + 208.*y.^3 + 96.*y.^2 - 144.*y + 32).*x.^3 + (288.*y.^4 - 312.*y.^3 - 144.*y.^2 + 216.*y - 48).*x.^2 + (- 96.*y.^4 + 104.*y.^3 + 48.*y.^2 - 72.*y + 16).*x) + x.*(- 96.*y.^3 + 96.*y.^2 + 8.*y - 16) + z.^2.*(x.^2.*(- 288.*y.^4 + 312.*y.^3 + 144.*y.^2 - 216.*y + 48) - 20.*y - x.^3.*(- 384.*y.^4 + 704.*y.^3 - 384.*y.^2 + 64.*y) + x.*(96.*y.^3 - 96.*y.^2 + 40.*y - 16) + 48.*y.^2 - 48.*y.^3 + 8) - z.^5.*(((432.*y.^3)/5 - (864.*y.^2)/5 + (528.*y)/5 - 96/5).*x.^3 + (- (648.*y.^3)/5 + (1296.*y.^2)/5 - (792.*y)/5 + 144/5).*x.^2 + ((216.*y.^3)/5 - (432.*y.^2)/5 + (264.*y)/5 - 48/5).*x))./(-w).^3) .* (y>0.5);
%! d2Fduv(3,:,:,:) = (((x.^2.*(48.*z.^2 - 48.*z.^4) - x.^4.*(16.*z.^4 - 48.*z.^2 + 32) + x.^3.*(48.*z.^4 - 96.*z.^2 + 32) + 16.*x.*z.^4).*y.^4 + (x.^2.*(- 48.*z.^5 + 48.*z.^3 + 144.*z.^2 - 144) - x.^3.*(- 32.*z.^5 + 32.*z.^3 + 48.*z.^2 - 48) + x.*(16.*z.^5 - 144.*z.^2 + 96) + 48.*z.^2).*y.^3 + (x.*(96.*z.^2 - 48) + x.^3.*(48.*z.^2 - 48) - x.^2.*(120.*z.^2 - 96) - 24.*z.^2).*y.^2 + (x.*(16.*z.^3 - 24) - 8.*z.^3 + 24).*y + 8.*x - 8)./w.^3) .* (y<0.5) + ...
%!     ((8.*y - x.^4.*(96.*y.^4 - 320.*y.^3 + 384.*y.^2 - 192.*y + 32) + x.^3.*(96.*y.^4 - 368.*y.^3 + 528.*y.^2 - 336.*y + 80) + z.^3.*((288.*y.^3 - 576.*y.^2 + 352.*y - 64).*x.^3 + (- 432.*y.^3 + 864.*y.^2 - 528.*y + 96).*x.^2 + (48.*y - 32).*x - 24.*y + 16) - x.*(96.*y.^3 - 240.*y.^2 + 200.*y - 56) - z.^4.*((48.*y.^4 - 160.*y.^3 + 192.*y.^2 - 96.*y + 16).*x.^4 + (- 144.*y.^4 + 480.*y.^3 - 576.*y.^2 + 288.*y - 48).*x.^3 + (144.*y.^4 - 480.*y.^3 + 576.*y.^2 - 288.*y + 48).*x.^2 + (- 48.*y.^4 + 160.*y.^3 - 192.*y.^2 + 96.*y - 16).*x) + z.^2.*(x.^4.*(144.*y.^4 - 480.*y.^3 + 576.*y.^2 - 288.*y + 48) - 96.*y + x.^2.*(144.*y.^4 - 624.*y.^3 + 984.*y.^2 - 672.*y + 168) - x.^3.*(288.*y.^4 - 1008.*y.^3 + 1296.*y.^2 - 720.*y + 144) + x.*(144.*y.^3 - 384.*y.^2 + 336.*y - 96) + 120.*y.^2 - 48.*y.^3 + 24) - z.^5.*((288.*y.^3 - 576.*y.^2 + 352.*y - 64).*x.^3 + (- 432.*y.^3 + 864.*y.^2 - 528.*y + 96).*x.^2 + (144.*y.^3 - 288.*y.^2 + 176.*y - 32).*x) + x.^2.*(144.*y.^3 - 384.*y.^2 + 336.*y - 96) - 8)./(-w).^3) .* (y>0.5);
%! d2Fduw = zeros ([3, size(x)]);
%! d2Fduw(1,:,:,:) = ((x.^2.*((24.*y.^4.*z.^2)/5 + 2.*z.*(8.*y.^4 + 4.*y.^2)) - (12.*y.^2.*z.^2)/5 + (24.*x.*y.^2.*z.^2)/5)./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(x.*((8.*y.^2.*z.^3)/5 + 8.*y.^2) - (4.*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8.*y.^4 + 4.*y.^2) + (8.*y.^4.*z.^3)/5 - 4.*y.^2) + 2))./w.^3) .* (y<0.5) + ...
%!     (-((- (4.*(3.*y - 1).*(y - 1).*(36.*y.^4 - 96.*y.^3 + 88.*y.^2 - 32.*y + 4).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(36.*y.^4 - 96.*y.^3 + 100.*y.^2 - 48.*y + 8).*x.^3)/5 - (4.*(3.*y - 1).*(y - 1).*(18.*y.^2 - 24.*y + 6).*x.^2)/5 + (4.*(3.*y - 1).*(y - 1).*(6.*y.^2 - 8.*y + 2).*x)/5).*z.^4 + ((4.*x.^3.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 820.*y.^2 - 240.*y + 20))/5 - (4.*x.^4.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 820.*y.^2 - 240.*y + 20))/5).*z.^3 + (- (4.*(3.*y - 1).*(y - 1).*(108.*y.^4 - 288.*y.^3 + 264.*y.^2 - 96.*y + 12).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(36.*y.^2 - 48.*y + 12).*x.^3)/5 - (24.*(3.*y - 1).*(y - 1).*x)/5 + (12.*(3.*y - 1).*(y - 1))/5).*z.^2 + (- (4.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 940.*y.^2 - 400.*y + 60).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(360.*y.^2 - 480.*y + 120).*x.^3)/5 - (4.*(3.*y - 1).*(y - 1).*(180.*y.^2 - 240.*y + 90).*x.^2)/5 + 16.*(3.*y - 1).*(y - 1).*x).*z)./(-w).^3) .* (y>0.5);
%! d2Fduw(2,:,:,:) = ((2.*z.*(8.*x.^2.*y.^4 - y.^2.*(8.*y - 4.*y.^2) + (2.*x.*y.^2.*(40.*y - 20.*y.^2))/5) + 3.*z.^2.*((12.*x.^2.*y.^4)/5 + (12.*x.*y.^2)/5 - (6.*y.^2)/5))./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(z.^2.*(8.*x.^2.*y.^4 - y.^2.*(8.*y - 4.*y.^2) + (2.*x.*y.^2.*(40.*y - 20.*y.^2))/5) + z.^3.*((12.*x.^2.*y.^4)/5 + (12.*x.*y.^2)/5 - (6.*y.^2)/5) + (2.*x.*y.^2.*(20.*y.^2 - 40.*y + 20))/5))./w.^3) .* (y<0.5) + ...
%!     (((6.*(3.*y.^2 - 4.*y + 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 6.*x + 3).*z.^2)/5 + (4.*(3.*y.^2 - 4.*y + 1).*(60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2 - 20.*x.*y.^2 - 10.*x + 10.*y.^2 + 5).*z)/5)./w.^2 - (2.*((2.*(3.*y.^2 - 4.*y + 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 6.*x + 3).*z.^3)/5 + (2.*(3.*y.^2 - 4.*y + 1).*(60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2 - 20.*x.*y.^2 - 10.*x + 10.*y.^2 + 5).*z.^2)/5 - (2.*(10.*x - 20.*x.*y.^2).*(3.*y.^2 - 4.*y + 1))/5).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fduw(3,:,:,:) = (- (2.*z.*(4.*y.*(2.*y.^2 - 3.*y.^3).*x.^2 - 4.*y.*(4.*y.^2 - 6.*y.^3).*x + 4.*y.*(2.*y.^2 - 3.*y.^3)) - 3.*z.^2.*(8.*x.^2.*y.^4 + 8.*x.*y.^2 - 4.*y.^2))./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(4.*y.*(3.*y - 2) + z.^3.*(8.*x.^2.*y.^4 + 8.*x.*y.^2 - 4.*y.^2) - z.^2.*(4.*y.*(2.*y.^2 - 3.*y.^3).*x.^2 - 4.*y.*(4.*y.^2 - 6.*y.^3).*x + 4.*y.*(2.*y.^2 - 3.*y.^3)) + 4.*x.^2.*y.*(4.*y.^2 - 6.*y.^3) - 4.*x.*y.*(- 6.*y.^3 + 4.*y.^2 + 3.*y - 2)))./w.^3) .* (y<0.5) + ...
%!     ((2.*z.*(4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1).*x.^2 - 4.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2).*x + 4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1)) + 3.*z.^2.*(4.*(y - 1).*(18.*y.^3 - 30.*y.^2 + 14.*y - 2).*x.^2 - 4.*(6.*y - 2).*(y - 1).*x + 4.*(3.*y - 1).*(y - 1)))./w.^2 - (2.*(z.^2.*(4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1).*x.^2 - 4.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2).*x + 4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1)) - 4.*(y - 1).^2 + z.^3.*(4.*(y - 1).*(18.*y.^3 - 30.*y.^2 + 14.*y - 2).*x.^2 - 4.*(6.*y - 2).*(y - 1).*x + 4.*(3.*y - 1).*(y - 1)) + 4.*x.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 11.*y - 3) - 4.*x.^2.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2)).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:,:) = (-(8.*x.*(x - 1).*(z.^3 + 5.*x.*z.^2 - 5.*x).*(- 6.*x.^2.*y.^2.*z.^2 + 6.*x.^2.*y.^2 + 6.*x.*y.^2.*z.^2 - 1))/5./w.^3) .* (y<0.5) + ...
%!     ((8.*x.*(x - 1).*(z.^3 + 5.*x.*z.^2 - 5.*x).*(- 54.*x.^2.*y.^2.*z.^2 + 54.*x.^2.*y.^2 + 72.*x.^2.*y.*z.^2 - 72.*x.^2.*y - 26.*x.^2.*z.^2 + 26.*x.^2 + 54.*x.*y.^2.*z.^2 - 72.*x.*y.*z.^2 + 26.*x.*z.^2 + 3))/5./(-w).^3) .* (y>0.5);
%! d2Fdvv(2,:,:,:) = ((2.*((8.*x.*z.^2 - x.^2.*(8.*z.^2 - 8)).*y.^2 + ((12.*x.*z.^3)/5 - x.^2.*((12.*z.^3)/5 + 8) + 4).*y - 4).*(- 4.*y.*x.^2.*z.^2 + 4.*y.*x.^2 + 4.*y.*x.*z.^2))./w.^3 - ((12.*x.*z.^3)/5 + 2.*y.*(8.*x.*z.^2 - x.^2.*(8.*z.^2 - 8)) - x.^2.*((12.*z.^3)/5 + 8) + 4)./w.^2) .* (y<0.5) + ...
%!     ((z.^2.*(x.*(32.*y + 4) - x.^2.*(32.*y + 4)) + x.^2.*(32.*y - 20) + z.^3.*((36.*x)/5 - (36.*x.^2)/5) + 4)./w.^2 - (2.*(4.*y + z.^3.*(x.*((36.*y)/5 - 24/5) - x.^2.*((36.*y)/5 - 24/5)) + z.^2.*(x.*(16.*y.^2 + 4.*y - 8) - x.^2.*(16.*y.^2 + 4.*y - 8)) + x.^2.*(16.*y.^2 - 20.*y + 8)).*(8.*x.^2.*z.^2 + 12.*x.^2.*y - 8.*x.*z.^2 - 8.*x.^2 + 12.*x.*y.*z.^2 - 12.*x.^2.*y.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fdvv(3,:,:,:) = ((2.*y.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(3.*x - 3).*(x - 1) + 8.*x.*z.^3.*(x - 1))./w.^2 - (2.*(4.*(x - 1).^2 - y.*(4.*(3.*x - 3).*(x - 1) - 8.*x.*z.^3.*(x - 1)) + y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1))).*(- 4.*y.*x.^2.*z.^2 + 4.*y.*x.^2 + 4.*y.*x.*z.^2))./w.^3) .* (y<0.5) + ...
%!     ((4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1) + 2.*y.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 24.*x.*z.^3.*(x - 1) - 4.*z.^2.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x))./w.^2 - (2.*(y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(x - 1).*(2.*x.^3 - 2.*x.^2 + x - 1) - y.*(24.*x.*(x - 1).*z.^3 + 4.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z.^2 - 4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1)) + 16.*x.*z.^3.*(x - 1) + 4.*z.^2.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x)).*(8.*x.^2.*z.^2 + 12.*x.^2.*y - 8.*x.*z.^2 - 8.*x.^2 + 12.*x.*y.*z.^2 - 12.*x.^2.*y.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fdvw = zeros ([3, size(x)]);
%! d2Fdvw(1,:,:,:) = (((8.*x.*z.*(x - 1).*(20.*x.^3.*z.^2 - 20.*x.^3 + 2.*x.^2.*z.^3 - 20.*x.^2.*z.^2 + 6.*x.^2.*z + 40.*x.^2 - 2.*x.*z.^3).*y.^3)/5 + (8.*x.*z.*(10.*x + 3.*z).*(x - 1).*y)/5)./w.^3) .* (y<0.5) + ...
%!     (((8.*x.*(3.*y - 2).*(x - 1).*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*z.^4)/5 + (8.*x.*(3.*y - 2).*(x - 1).*(- 60.*x.^3.*y.^2 + 80.*x.^3.*y - 20.*x.^3 + 60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2).*z.^3)/5 - (8.*x.*(3.*y - 2).*(x - 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 3).*z.^2)/5 + (8.*x.*(3.*y - 2).*(x - 1).*(60.*x.^3.*y.^2 - 80.*x.^3.*y + 20.*x.^3 - 120.*x.^2.*y.^2 + 160.*x.^2.*y - 40.*x.^2 + 10.*x).*z)/5)./(-w).^3) .* (y>0.5);
%! d2Fdvw(2,:,:,:) = ((4.*x.*y.*z.*(x - 1).*(40.*x.^2.*y.^3.*z.^2 - 40.*x.^2.*y.^3 + 6.*x.^2.*y.^2.*z.^3 + 18.*x.^2.*y.^2.*z + 80.*x.^2.*y.^2 - 40.*x.*y.^3.*z.^2 - 6.*x.*y.^2.*z.^3 - 40.*y.^2 + 60.*y + 9.*z))/5./w.^3) .* (y<0.5) + ...
%!     (-((4.*x.*(x - 1).*(54.*x.^2.*y.^3 - 108.*x.^2.*y.^2 + 66.*x.^2.*y - 12.*x.^2 - 54.*x.*y.^3 + 108.*x.*y.^2 - 66.*x.*y + 12.*x).*z.^4)/5 + (4.*x.*(x - 1).*(240.*x.^2.*y.^4 - 260.*x.^2.*y.^3 - 120.*x.^2.*y.^2 + 180.*x.^2.*y - 40.*x.^2 - 240.*x.*y.^4 + 260.*x.*y.^3 + 120.*x.*y.^2 - 180.*x.*y + 40.*x).*z.^3)/5 - (4.*x.*(x - 1).*(- 162.*x.^2.*y.^3 + 324.*x.^2.*y.^2 - 198.*x.^2.*y + 36.*x.^2 + 27.*y - 18).*z.^2)/5 - (4.*x.*(x - 1).*(240.*x.^2.*y.^4 - 980.*x.^2.*y.^3 + 1320.*x.^2.*y.^2 - 700.*x.^2.*y + 120.*x.^2 + 120.*y.^3 - 120.*y.^2 + 50.*y - 20).*z)/5)./(-w).^3) .* (y>0.5);
%! d2Fdvw(3,:,:,:) = (-(y.^3.*(8.*x.*z.*(x - 1).*(12.*x.^2 - 24.*x + 12) - 48.*x.^3.*z.^2.*(x - 1) + 8.*x.*z.^4.*(2.*x - 2.*x.^2).*(x - 1)) + y.^4.*(8.*x.*(x - 1).*(- 4.*x.^4 + 12.*x.^3 - 12.*x.^2 + 4.*x).*z.^3 + 8.*x.*(x - 1).*(4.*x.^4 - 8.*x.^3 + 4.*x.^2).*z) - 24.*x.*y.*z.^2.*(x - 1) - 8.*x.*y.^2.*z.*(x - 1).*(6.*x.^2 - 12.*x + 6))./w.^3) .* (y<0.5) + ...
%!     ((8.*z.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x) - y.*(72.*x.*(x - 1).*z.^2 + 8.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z) + 48.*x.*z.^2.*(x - 1) + 8.*y.^2.*z.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x))./w.^2 - (2.*(y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(x - 1).*(2.*x.^3 - 2.*x.^2 + x - 1) - y.*(24.*x.*(x - 1).*z.^3 + 4.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z.^2 - 4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1)) + 16.*x.*z.^3.*(x - 1) + 4.*z.^2.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x)).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fdww = zeros ([3, size(x)]);
%! d2Fdww(1,:,:,:) = ((32.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(5.*x + z + 10.*x.^2.*y.^2 + 2.*x.^2.*y.^2.*z))./(5.*w.^3) - (8.*x.*y.^2.*(x - 1).*(15.*x + z + 30.*x.^2.*y.^2 + 2.*x.^2.*y.^2.*z))/5./w.^2) .* (y<0.5) + ...
%!     (((8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(36.*x.^4.*y.^4 - 96.*x.^4.*y.^3 + 88.*x.^4.*y.^2 - 32.*x.^4.*y + 4.*x.^4 - 36.*x.^3.*y.^4 + 96.*x.^3.*y.^3 - 88.*x.^3.*y.^2 + 32.*x.^3.*y - 4.*x.^3 - 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*z.^3)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(540.*x.^4.*y.^4 - 1440.*x.^4.*y.^3 + 1320.*x.^4.*y.^2 - 480.*x.^4.*y + 60.*x.^4 - 540.*x.^3.*y.^4 + 1440.*x.^3.*y.^3 - 1410.*x.^3.*y.^2 + 600.*x.^3.*y - 90.*x.^3 + 90.*x.^2.*y.^2 - 120.*x.^2.*y + 30.*x.^2).*z.^2)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(108.*x.^4.*y.^4 - 288.*x.^4.*y.^3 + 264.*x.^4.*y.^2 - 96.*x.^4.*y + 12.*x.^4 - 36.*x.^2.*y.^2 + 48.*x.^2.*y - 12.*x.^2 + 3).*z)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(180.*x.^4.*y.^4 - 480.*x.^4.*y.^3 + 440.*x.^4.*y.^2 - 160.*x.^4.*y + 20.*x.^4 - 30.*x.^3.*y.^2 + 40.*x.^3.*y - 10.*x.^3 - 30.*x.^2.*y.^2 + 40.*x.^2.*y - 10.*x.^2 + 5.*x))/5)./(-w).^3) .* (y>0.5);
%! d2Fdww(2,:,:,:) = ((16.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(20.*y + 3.*z + 20.*x.^2.*y.^2 - 10.*y.^2 + 6.*x.^2.*y.^2.*z))./(5.*w.^3) - (12.*x.*y.^2.*(x - 1).*(20.*y + z + 20.*x.^2.*y.^2 - 10.*y.^2 + 2.*x.^2.*y.^2.*z))/5./w.^2) .* (y<0.5) + ...
%!     (((4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(108.*x.^4.*y.^4 - 288.*x.^4.*y.^3 + 264.*x.^4.*y.^2 - 96.*x.^4.*y + 12.*x.^4 - 108.*x.^3.*y.^4 + 288.*x.^3.*y.^3 - 264.*x.^3.*y.^2 + 96.*x.^3.*y - 12.*x.^3 - 18.*x.^2.*y.^2 + 24.*x.^2.*y - 6.*x.^2 + 18.*x.*y.^2 - 24.*x.*y + 6.*x).*z.^3)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(1080.*x.^4.*y.^4 - 2880.*x.^4.*y.^3 + 2640.*x.^4.*y.^2 - 960.*x.^4.*y + 120.*x.^4 - 1080.*x.^3.*y.^4 + 2880.*x.^3.*y.^3 - 2640.*x.^3.*y.^2 + 960.*x.^3.*y - 120.*x.^3 - 180.*x.^2.*y.^4 + 240.*x.^2.*y.^3 - 150.*x.^2.*y.^2 + 120.*x.^2.*y - 30.*x.^2 + 180.*x.*y.^4 - 240.*x.*y.^3 + 150.*x.*y.^2 - 120.*x.*y + 30.*x).*z.^2)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(324.*x.^4.*y.^4 - 864.*x.^4.*y.^3 + 792.*x.^4.*y.^2 - 288.*x.^4.*y + 36.*x.^4 - 108.*x.^2.*y.^2 + 144.*x.^2.*y - 36.*x.^2 + 9).*z)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(360.*x.^4.*y.^4 - 960.*x.^4.*y.^3 + 880.*x.^4.*y.^2 - 320.*x.^4.*y + 40.*x.^4 - 60.*x.^2.*y.^4 + 80.*x.^2.*y.^3 - 110.*x.^2.*y.^2 + 120.*x.^2.*y - 30.*x.^2 + 10.*y.^2 + 5))/5)./(-w).^3) .* (y>0.5);
%! d2Fdww(3,:,:,:) = ((32.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(2.*y + z - 3.*x.^2.*y.^2 - 4.*x.*y + 6.*x.*y.^2 + 2.*x.^2.*y - 3.*y.^2 + 2.*x.^2.*y.^2.*z))./w.^3 - (8.*x.*y.^2.*(x - 1).*(6.*y + z - 9.*x.^2.*y.^2 - 12.*x.*y + 18.*x.*y.^2 + 6.*x.^2.*y - 9.*y.^2 + 2.*x.^2.*y.^2.*z))./w.^2) .* (y<0.5) + ...
%!    ((2.*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*(2.*z - 4.*y - 4.*x + 2.*x.^2.*y.^2 + 8.*x.*y - 4.*x.*y.^2 - 4.*x.^2.*y - 4.*x.^2.*z + 2.*x.^2 + 2.*y.^2 + 16.*x.^2.*y.*z - 12.*x.^2.*y.^2.*z + 2))./w.^2 - (8.*z.^2.*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).^2.*(2.*z - 4.*y - 4.*x + 2.*x.^2.*y.^2 + 8.*x.*y - 4.*x.*y.^2 - 4.*x.^2.*y - 4.*x.^2.*z + 2.*x.^2 + 2.*y.^2 + 16.*x.^2.*y.*z - 12.*x.^2.*y.^2.*z + 2))./(-w).^3 - (4.*z.*(12.*x.^2.*y.^2 - 16.*x.^2.*y + 4.*x.^2 - 2).*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x))./w.^2) .* (y>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (dFdu, jac{1}, 1e3*eps)
%! assert (dFdv, jac{2}, 1e3*eps)
%! assert (dFdw, jac{3}, 1e3*eps)
%! assert (d2Fduu, hess{1,1}, 1e3*eps)
%! assert (d2Fduv, hess{1,2}, 1e3*eps)
%! assert (d2Fduw, hess{1,3}, 1e3*eps)
%! assert (d2Fduv, hess{2,1}, 1e3*eps)
%! assert (d2Fdvv, hess{2,2}, 1e3*eps)
%! assert (d2Fdvw, hess{2,3}, 1e3*eps)
%! assert (d2Fduw, hess{3,1}, 1e3*eps)
%! assert (d2Fdvw, hess{3,2}, 1e3*eps)
%! assert (d2Fdww, hess{3,3}, 1e3*eps)



%!test
%! nrb = nrbextrude (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [0 0 1]);
%! nrb = nrbdegelev (nrb, [1 1 1]);
%! nrb.coefs (4,2,2,2) = 1.1;
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 24); Y = linspace (0, 1, 24); Z = linspace (0, 1, 24);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y Z});
%! [y, x, z] = meshgrid (X, Y, Z);
%! F = zeros ([3, size(x)]);
%! F(1,:,:,:) = (5.*x)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! F(2,:,:,:) = (5.*y)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! F(3,:,:,:) = (5.*z)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! dFdu = zeros ([3, size(x)]);
%! dFdu(1,:,:,:) = ((z.*(20.*y - 20.*y.^2) - z.^2.*(20.*y - 20.*y.^2)).*x.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! dFdu(2,:,:,:) = (y.^2.*(5.*z.*(8.*x - 4) - 5.*z.^2.*(8.*x - 4)) - y.^3.*(5.*z.*(8.*x - 4) - 5.*z.^2.*(8.*x - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdu(3,:,:,:) = (z.^2.*(5.*y.*(8.*x - 4) - 5.*y.^2.*(8.*x - 4)) - z.^3.*(5.*y.*(8.*x - 4) - 5.*y.^2.*(8.*x - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdv = zeros ([3, size(x)]);
%! dFdv(1,:,:,:) = (x.^2.*(5.*z.*(8.*y - 4) - 5.*z.^2.*(8.*y - 4)) - x.^3.*(5.*z.*(8.*y - 4) - 5.*z.^2.*(8.*y - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdv(2,:,:,:) = ((z.*(20.*x - 20.*x.^2) - z.^2.*(20.*x - 20.*x.^2)).*y.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! dFdv(3,:,:,:) = (z.^2.*(5.*x.*(8.*y - 4) - 5.*x.^2.*(8.*y - 4)) - z.^3.*(5.*x.*(8.*y - 4) - 5.*x.^2.*(8.*y - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw = zeros ([3, size(x)]);
%! dFdw(1,:,:,:) = (x.^2.*(y.*(40.*z - 20) - y.^2.*(40.*z - 20)) - x.^3.*(y.*(40.*z - 20) - y.^2.*(40.*z - 20)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw(2,:,:,:) = (y.^2.*(x.*(40.*z - 20) - x.^2.*(40.*z - 20)) - y.^3.*(x.*(40.*z - 20) - x.^2.*(40.*z - 20)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw(3,:,:,:) = ((y.*(20.*x - 20.*x.^2) - y.^2.*(20.*x - 20.*x.^2)).*z.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:,:) = (40.*y.*z.*(y - 1).*(z - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z + 15.*x - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduu(2,:,:,:) = (40.*y.^2.*z.*(y - 1).*(z - 1).*(4.*y.^2.*z.^2 - 4.*y.^2.*z - 4.*y.*z.^2 + 4.*y.*z + 5) - 40.*x.*y.^2.*z.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z) + 40.*x.^2.*y.^2.*z.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduu(3,:,:,:) = (40.*y.*z.^2.*(y - 1).*(z - 1).*(4.*y.^2.*z.^2 - 4.*y.^2.*z - 4.*y.*z.^2 + 4.*y.*z + 5) - 40.*x.*y.*z.^2.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z) + 40.*x.^2.*y.*z.^2.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:,:) = (20.*x.*z.*(2.*y - 1).*(z - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 15.*x - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv(2,:,:,:) = (20.*y.*z.*(2.*x - 1).*(z - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z + 15.*y - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv(3,:,:,:) = (20.*z.^2.*(2.*x - 1).*(2.*y - 1).*(z - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw = zeros ([3, size(x)]);
%! d2Fduw(1,:,:,:) = (20.*x.*y.*(2.*z - 1).*(y - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 15.*x - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw(2,:,:,:) = (20.*y.^2.*(2.*x - 1).*(2.*z - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw(3,:,:,:) = (20.*y.*z.*(2.*x - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.*z.^3 + 4.*x.^2.*y.*z.^2 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.^2.*z.^2 + 4.*x.*y.*z.^3 - 4.*x.*y.*z.^2 + 15.*z - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:,:) = (40.*x.^2.*z.*(x - 1).*(z - 1).*(4.*x.^2.*z.^2 - 4.*x.^2.*z - 4.*x.*z.^2 + 4.*x.*z + 5) + 40.*x.^2.*y.^2.*z.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z) - 40.*x.^2.*y.*z.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv(2,:,:,:) = (40.*x.*z.*(x - 1).*(z - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 15.*y - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv(3,:,:,:) = (40.*x.*z.^2.*(x - 1).*(z - 1).*(4.*x.^2.*z.^2 - 4.*x.^2.*z - 4.*x.*z.^2 + 4.*x.*z + 5) + 40.*x.*y.^2.*z.^2.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z) - 40.*x.*y.*z.^2.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw = zeros ([3, size(x)]);
%! d2Fdvw(1,:,:,:) = (20.*x.^2.*(2.*y - 1).*(2.*z - 1).*(x - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw(2,:,:,:) = (20.*x.*y.*(2.*z - 1).*(x - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z + 15.*y - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw(3,:,:,:) = (20.*x.*z.*(2.*y - 1).*(x - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.*z.^3 + 4.*x.^2.*y.*z.^2 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.^2.*z.^2 + 4.*x.*y.*z.^3 - 4.*x.*y.*z.^2 + 15.*z - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww = zeros ([3, size(x)]);
%! d2Fdww(1,:,:,:) = (40.*x.^2.*y.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y + 5) + 40.*x.^2.*y.*z.^2.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y) - 40.*x.^2.*y.*z.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww(2,:,:,:) = (40.*x.*y.^2.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y + 5) + 40.*x.*y.^2.*z.^2.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y) - 40.*x.*y.^2.*z.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww(3,:,:,:) = (40.*x.*y.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.*z.^3 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.*z.^3 + 15.*z - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;


%!test
%! knots = [0 0 0 0.5 1 1 1];
%! coefs(1,:) = [0 2 4 2];
%! coefs(2,:) = [0 2 2 0];
%! coefs(3,:) = [0 4 2 0];
%! coefs(4,:) = [1 2 2 1];
%! nrb = nrbmak (coefs, knots);
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! x = linspace (0, 1, 10);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, x);
%! w = -4*x.^2 + 4*x + 1;
%! F = zeros (3,numel(x)); DF = zeros (3, numel(x)); D2F = zeros (3, numel(x));
%! F(1,:) = (-4*x.*(x-2)./w) .* (x<0.5) + ((4*x - 5)./w + 3) .* (x>0.5);
%! F(2,:) = (2-2./w);
%! F(3,:) = (-4*x.*(5*x-4)./w) .* (x<0.5) + (-4*(x.^2 - 1)./w) .* (x>0.5);
%! DF(1,:) = (8*(2*x.^2-x+1)./w.^2) .* (x<0.5) + (8*(2*x-3).*(x-1)./w.^2) .* (x>0.5);
%! DF(2,:) = -8*(2*x-1)./w.^2;
%! DF(3,:) = -(8*(2*x.^2+5*x-2)./w.^2) .* (x<0.5) - (8*(2*x.^2-3*x+2)./w.^2) .* (x>0.5);
%! D2F(1,:) = 8*(16*x.^3-12*x.^2+24*x-9)./w.^3 .* (x<0.5) + 8*(16*x.^3-60*x.^2+72*x-29)./w.^3 .* (x>0.5);
%! D2F(2,:) = -16*(12*x.^2-12*x+5)./w.^3;
%! D2F(3,:) = -8*(16*x.^3+60*x.^2-48*x+21)./w.^3 .* (x<0.5) -8*(16*x.^3-36*x.^2+48*x-19)./w.^3 .* (x>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (DF, jac, 1e3*eps)
%! assert (D2F, hess, 1e3*eps)

%!test
%! knots = {[0 0 0 1 1 1], [0 0 0 0.5 1 1 1]};
%! coefs = ones (4,3,4);
%! coefs(1,:,:) = reshape ([0 0 0 0; 1 1 1 1; 2 2 4 2], 1, 3, 4);
%! coefs(2,:,:) = reshape ([0 1 2 3; 0 1 2 3; 0 1 4 3], 1, 3, 4);
%! coefs(3,:,:) = reshape ([0 1 0 0; 0 0 0 0; 0 0 0 0], 1, 3, 4);
%! coefs(4,:,:) = reshape ([1 1 1 1; 1 1 1 1; 1 1 2 1], 1, 3, 4);
%! nrb = nrbmak (coefs, knots);
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 4); Y = linspace (0, 1, 4);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y});
%! [y, x] = meshgrid (X, Y);
%! w = (2*x.^2.*y.^2 + 1) .* (y < 0.5) + (-6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 1) .* (y > 0.5);
%! F = zeros ([3,size(x)]);
%! F(1,:,:) = ((2*x - 2) ./w + 2) .* (y<0.5) + (2 + (2*x-2)./w) .* (y > 0.5);
%! F(2,:,:) = (2 - (2*(y-1).^2)./w).*(y<0.5) + ...
%!     ((-12*x.^2.*y.^2 + 16*x.^2.*y - 4*x.^2 + 2*y.^2 + 1)./w).*(y>0.5);
%! F(3,:,:) = (-2*y.*(3*y - 2).*(x - 1).^2./w) .* (y<0.5) + ...
%!     (2*(x - 1).^2.*(y - 1).^2./w) .* (y>0.5); 
%! dFdu = zeros ([3,size(x)]);
%! dFdu(1,:,:) = (((8*x - 4*x.^2).*y.^2 + 2)./w.^2).*(y<0.5) + ...
%!     (((12*y.^2 - 16*y + 4).*x.^2 + (-24*y.^2 + 32*y - 8).*x + 2)./w.^2).*(y>0.5);
%! dFdu(2,:,:) = (8*x.*y.^2.*(y - 1).^2./w.^2).*(y<0.5) + ...
%!     ((4*x.*(3*y - 1).*(2*y.^2 - 1).*(y - 1))./w.^2).*(y>0.5);
%! dFdu(3,:,:) = (-4*y.*(2.*x.*y.^2 + 1).*(3*y - 2).*(x - 1)./w.^2).*(y<0.5) + ...
%!     ((-4*(x - 1).*(y - 1).^2.*(6*x.*y.^2 - 8*x.*y + 2*x - 1))./w.^2).*(y>0.5);
%! dFdv = zeros ([3,size(x)]);
%! dFdv(1,:,:) = (-8*x.^2.*y.*(x - 1)./w.^2).*(y<0.5) + ...
%!     (8*x.^2.*(3*y - 2).*(x - 1)./w.^2).*(y>0.5);
%! dFdv(2,:,:) = (-4*(2*y.*x.^2 + 1).*(y - 1)./w.^2).*(y<0.5) + ...
%!     (((16*y.^2 - 20*y + 8).*x.^2 + 4*y)./w.^2).*(y>0.5);
%! dFdv(3,:,:) = (-4*(x - 1).^2.*(2*x.^2.*y.^2 + 3*y - 1)./w.^2).*(y<0.5) + ...
%!     (4*(x - 1).^2.*(y - 1).*(2*x.^2 - 2*x.^2.*y + 1)./w.^2).*(y>0.5);
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:) = (-((48*x.^2 - 16*x.^3).*y.^4 + (24*x - 8).*y.^2)./w.^3).*(y<0.5) + ...
%!     (((32*(3*y - 1).*(x - 1).*(y - 1))-(8*(3*y - 1).*(x - 3).*(y - 1).*w))./w.^3).*(y>0.5);
%! d2Fduu(2,:,:) = (-(8*y.^2.*(6*x.^2.*y.^2 - 1).*(y - 1).^2)./w.^3).*(y<0.5) + ...
%!     ((4*(3*y - 1).*(2*y.^2 - 1).*(y - 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 + 1))./w.^3).*(y>0.5);
%! d2Fduu(3,:,:) =  ((4*y.*(3*y - 2).*(8*x.^3.*y.^4 - 12*x.^2.*y.^4 + 6*x.^2.*y.^2 - 12*x.*y.^2 + 2*y.^2 - 1))./w.^3).*(y<0.5) + ...
%!  ((4*(y - 1).^2.*(6*y.^2 - 8*y + 3) - 4*x.^3.*(y - 1).^2.*(72*y.^4 - 192*y.^3 + 176*y.^2 - 64*y + 8) + 4*x.^2.*(y - 1).^2.*(108*y.^4 - 288*y.^3 + 282*y.^2 - 120*y + 18) - 4*x.*(y - 1).^2.*(36*y.^2 - 48*y + 12))./w.^3) .* (y>0.5);
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:) = (8*x.^2.*(6*x.^2.*y.^2 - 1).*(x - 1)./w.^3) .* (y<0.5) + ...
%!      (8*x.^2.*(x - 1).*(54*x.^2.*y.^2 - 72*x.^2.*y + 26*x.^2 + 3)./w.^3) .* (y>0.5);
%! d2Fdvv(2,:,:) =  (-((48*y.^2 - 32*y.^3).*x.^4 + (- 24*y.^2 + 48*y - 8).*x.^2 + 4)./w.^3) .*(y<0.5) + ...
%!       (((192*y.^3 - 360*y.^2 + 288*y - 88).*x.^4 + (72*y.^2 - 28).*x.^2 + 4)./w.^3) .* (y>0.5);
%! d2Fdvv(3,:,:) =  (4*(x - 1).^2.*(8*x.^4.*y.^3 + 18*x.^2.*y.^2 - 12*x.^2.*y - 3))./w.^3 .* (y<0.5) + ...
%!     ((4*(x - 1).^2.*(24*x.^4 + 18*x.^2 + 1) + 4*y.^2.*(72*x.^4 + 18*x.^2).*(x - 1).^2 - 96*x.^4.*y.^3.*(x - 1).^2 - 4*y.*(72*x.^4 + 36*x.^2).*(x - 1).^2)./w.^3) .* (y>0.5);
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:) = (-(y.^3.*(32*x.^3 - 16*x.^4) - y.*(16*x - 24*x.^2))./w.^3) .* (y<0.5) + ...
%!     (-(-8*(3*y - 2).*(6*y.^2 - 8*y + 2).*x.^4 + 8*(3*y - 2).*(12*y.^2 - 16*y + 4).*x.^3 + (48 - 72*y).*x.^2 + (48*y - 32).*x)./w.^3) .* (y>0.5);
%! d2Fduv(2,:,:) = (16*x.*y.*(y - 1).*(2*x.^2.*y.^2 + 2*y - 1)./w.^3) .* (y<0.5) + ...
%!     (-(8*x.*(4*y.^2 - 5*y + 2))./w.^2 + (16*x.*(3*y - 2).*(2*y.^2 - 1))./w.^3) .* (y>0.5);
%! d2Fduv(3,:,:) = (-(8*(x - 1).*(4*x.^3.*y.^4 - 6*x.^2.*y.^3 + 6*x.^2.*y.^2 + 12*x.*y.^3 - 6*x.*y.^2 + 3*y - 1))./w.^3) .* (y<0.5) + ...
%!     ((8*(x - 1).*(y - 1).*(12*x.^3.*y.^3 - 28*x.^3.*y.^2 + 20*x.^3.*y - 4*x.^3 + 6*x.^2.*y.^2 - 12*x.^2.*y + 6*x.^2 - 12*x.*y.^2 + 18*x.*y - 6*x + 1))./w.^3) .* (y>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (dFdu, jac{1}, 1e3*eps)
%! assert (dFdv, jac{2}, 1e3*eps)
%! assert (d2Fduu, hess{1,1}, 1e3*eps)
%! assert (d2Fduv, hess{1,2}, 1e3*eps)
%! assert (d2Fduv, hess{2,1}, 1e3*eps)
%! assert (d2Fdvv, hess{2,2}, 1e3*eps)

%!test
%! knots = {[0 0 0 1 1 1], [0 0 0 0.5 1 1 1]};
%! coefs = ones (4,3,4);
%! coefs(1,:,:) = reshape ([0 0 0 0; 1 1 1 1; 2 2 4 2], 1, 3, 4);
%! coefs(2,:,:) = reshape ([0 1 2 3; 0 1 2 3; 0 1 4 3], 1, 3, 4);
%! coefs(3,:,:) = reshape ([0 1 0 0; 0 0 0 0; 0 0 0 0], 1, 3, 4);
%! coefs(4,:,:) = reshape ([1 1 1 1; 1 1 1 1; 1 1 2 1], 1, 3, 4);
%! nrb = nrbmak (coefs, knots);
%! nrb = nrbdegelev (nrbextrude (nrb, [0.4 0.6 2]), [0 0 1]);
%! nrb.coefs(4,2,3,3) = 1.5;
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 4); Y = linspace (0, 1, 4); Z = linspace (0, 1, 4);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y Z});
%! [y, x, z] = meshgrid (X, Y, Z);
%! w = (-2*x.^2.*y.^2.*z.^2 + 2*x.^2.*y.^2 + 2*x.*y.^2.*z.^2 + 1) .* (y < 0.5) + ...
%!     (6*x.^2.*y.^2.*z.^2 - 6*x.^2.*y.^2 - 8*x.^2.*y.*z.^2 + 8*x.^2.*y + 2*x.^2.*z.^2 - 2*x.^2 - 6*x.*y.^2.*z.^2 + 8*x.*y.*z.^2 - 2*x.*z.^2 + 1) .* (y > 0.5);
%! F = zeros ([3,size(x)]);
%! F(1,:,:,:) = ((10*x + 20*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2))./(5*w)) .* (y<0.5) + ...
%!     (60*x.^2.*y.^2 - 10*x + z.*(12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2) - 80*x.^2.*y + 20*x.^2)./(-5*w) .* (y > 0.5);
%! F(2,:,:,:) = ((20*y + 20*x.^2.*y.^2 + z.*(6*x.^2.*y.^2 + 3) - 10*y.^2)./(5*w)).*(y<0.5) + ...
%!     ((60*x.^2.*y.^2 + z.*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3) - 80*x.^2.*y + 20*x.^2 - 10*y.^2 - 5)./(-5*w)).*(y>0.5);
%! F(3,:,:,:) = ((4*y - 6*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2) - 8*x.*y + 12*x.*y.^2 + 4*x.^2.*y - 6*y.^2)./w) .* (y<0.5) + ...
%!     ((2*z - 4*y - 4*x + 2*x.^2.*y.^2 + 8*x.*y - 4*x.*y.^2 - 4*x.^2.*y - 4*x.^2.*z + 2*x.^2 + 2*y.^2 + 16*x.^2.*y.*z - 12*x.^2.*y.^2.*z + 2)./w) .* (y>0.5); 
%! dFdu = zeros ([3,size(x)]);
%! dFdu(1,:,:,:) = ((x.*((8*y.^2.*z.^3)/5 + 8*y.^2) - (4*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 2)./w.^2).*(y<0.5) + ...
%!     ((z.^3.*(x.^2.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (16*y)/5 - x.*((24*y.^2)/5 - (32*y)/5 + 8/5) + (12*y.^2)/5 + 4/5) - x.*(24*y.^2 - 32*y + 8) + x.^2.*(12*y.^2 - 16*y + 4) + x.^2.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) + 2)./w.^2).*(y>0.5);
%! dFdu(2,:,:,:) = ((z.^2.*(8*x.^2.*y.^4 - y.^2.*(8*y - 4*y.^2) + (2*x.*y.^2.*(40*y - 20*y.^2))/5) + z.^3.*((12*x.^2.*y.^4)/5 + (12*x.*y.^2)/5 - (6*y.^2)/5) + (2*x.*y.^2.*(20*y.^2 - 40*y + 20))/5)./w.^2).*(y<0.5) + ...
%!     (((2*(3*y.^2 - 4*y + 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 6*x + 3).*z.^3)/5 + (2*(3*y.^2 - 4*y + 1).*(60*x.^2.*y.^2 - 80*x.^2.*y + 20*x.^2 - 20*x.*y.^2 - 10*x + 10*y.^2 + 5).*z.^2)/5 - (2*(10*x - 20*x.*y.^2).*(3*y.^2 - 4*y + 1))/5)./w.^2).*(y>0.5);
%! dFdu(3,:,:,:) = ((4*y.*(3*y - 2) + z.^3.*(8*x.^2.*y.^4 + 8*x.*y.^2 - 4*y.^2) - z.^2.*(4*y.*(2*y.^2 - 3*y.^3).*x.^2 - 4*y.*(4*y.^2 - 6*y.^3).*x + 4*y.*(2*y.^2 - 3*y.^3)) + 4*x.^2.*y.*(4*y.^2 - 6*y.^3) - 4*x.*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2)) ./w.^2).*(y<0.5) + ...
%!     ((z.^2.*(4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1).*x.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2).*x + 4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)) - 4*(y - 1).^2 + z.^3.*(4*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2).*x.^2 - 4*(6*y - 2).*(y - 1).*x + 4*(3*y - 1).*(y - 1)) + 4*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) - 4*x.^2.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2))./w.^2) .* (y > 0.5); 
%! dFdv = zeros ([3,size(x)]);
%! dFdv(1,:,:,:) = ((8*x.*y.*(x - 1).*(z.^3 + 5*x.*z.^2 - 5*x))/5./w.^2).*(y<0.5) + ...
%!     (-(8*x.*(3*y - 2).*(x - 1).*(z.^3 + 5*x.*z.^2 - 5*x))/5./w.^2).*(y>0.5);
%! dFdv(2,:,:,:) = (-((8*x.*z.^2 - x.^2.*(8*z.^2 - 8)).*y.^2 + ((12*x.*z.^3)/5 - x.^2.*((12*z.^3)/5 + 8) + 4).*y - 4)./w.^2).*(y<0.5) + ...
%!     ((4*y + z.^3.*(x.*((36*y)/5 - 24/5) - x.^2.*((36*y)/5 - 24/5)) + z.^2.*(x.*(16*y.^2 + 4*y - 8) - x.^2.*(16*y.^2 + 4*y - 8)) + x.^2.*(16*y.^2 - 20*y + 8))./w.^2).*(y>0.5);
%! dFdv(3,:,:,:) = ((4*(x - 1).^2 - y.*(4*(3*x - 3).*(x - 1) - 8*x.*z.^3.*(x - 1)) + y.^2.*(4*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x).*z.^2 + 4*(2*x.^2 - 2*x.^3).*(x - 1)))./w.^2).*(y<0.5) + ...
%!     ((y.^2.*(4*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x).*z.^2 + 4*(2*x.^2 - 2*x.^3).*(x - 1)) - 4*(x - 1).*(2*x.^3 - 2*x.^2 + x - 1) - y.*(24*x.*(x - 1).*z.^3 + 4*(x - 1).*(4*x.^3 - 8*x.^2 + 4*x).*z.^2 - 4*(x - 1).*(4*x.^3 - 4*x.^2 + x - 1)) + 16*x.*z.^3.*(x - 1) + 4*z.^2.*(x - 1).*(2*x.^3 - 4*x.^2 + 2*x))./w.^2).*(y>0.5);
%! dFdw = zeros ([3,size(x)]);
%! dFdw(1,:,:,:) = ((4*x.^2.*y.^2 + 2)./(- 10*x.^2.*y.^2.*z.^2 + 10*x.^2.*y.^2 + 10*x.*y.^2.*z.^2 + 5) - ((20*x.*y.^2.*z - 20*x.^2.*y.^2.*z).*(10*x + 20*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2)))./(5*w).^2).*(y<0.5) + ...
%!     ((12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2)./(- 30*x.^2.*y.^2.*z.^2 + 30*x.^2.*y.^2 + 40*x.^2.*y.*z.^2 - 40*x.^2.*y - 10*x.^2.*z.^2 + 10*x.^2 + 30*x.*y.^2.*z.^2 - 40*x.*y.*z.^2 + 10*x.*z.^2 - 5) - ((60*x.^2.*y.^2 - 10*x + z.*(12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2) - 80*x.^2.*y + 20*x.^2).*(- 60*z.*x.^2.*y.^2 + 80*z.*x.^2.*y - 20*z.*x.^2 + 60*z.*x.*y.^2 - 80*z.*x.*y + 20*z.*x))./(5*w).^2).*(y>0.5);
%! dFdw(2,:,:,:) = ((6*x.^2.*y.^2 + 3)./(- 10*x.^2.*y.^2.*z.^2 + 10*x.^2.*y.^2 + 10*x.*y.^2.*z.^2 + 5) - ((20*x.*y.^2.*z - 20*x.^2.*y.^2.*z).*(20*y + 20*x.^2.*y.^2 + z.*(6*x.^2.*y.^2 + 3) - 10*y.^2))./(5*w).^2).*(y<0.5) + ...
%!     ((18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3)./(- 30*x.^2.*y.^2.*z.^2 + 30*x.^2.*y.^2 + 40*x.^2.*y.*z.^2 - 40*x.^2.*y - 10*x.^2.*z.^2 + 10*x.^2 + 30*x.*y.^2.*z.^2 - 40*x.*y.*z.^2 + 10*x.*z.^2 - 5) - ((- 60*z.*x.^2.*y.^2 + 80*z.*x.^2.*y - 20*z.*x.^2 + 60*z.*x.*y.^2 - 80*z.*x.*y + 20*z.*x).*(60*x.^2.*y.^2 + z.*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 3) - 80*x.^2.*y + 20*x.^2 - 10*y.^2 - 5))./(5*w).^2).*(y>0.5);
%! dFdw(3,:,:,:) = ((4*x.^2.*y.^2 + 2)./(2*x.^2.*y.^2 - z.^2.*(2*x.^2.*y.^2 - 2*x.*y.^2) + 1) + (2*z.*(2*x.^2.*y.^2 - 2*x.*y.^2).*(4*y - 6*x.^2.*y.^2 + z.*(4*x.^2.*y.^2 + 2) - 8*x.*y + 12*x.*y.^2 + 4*x.^2.*y - 6*y.^2))./w.^2).*(y<0.5) + ...
%!     ((12*x.^2.*y.^2 - 16*x.^2.*y + 4*x.^2 - 2)./(6*x.^2.*y.^2 + z.^2.*(- 6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 6*x.*y.^2 - 8*x.*y + 2*x) - 8*x.^2.*y + 2*x.^2 - 1) + (2*z.*(- 6*x.^2.*y.^2 + 8*x.^2.*y - 2*x.^2 + 6*x.*y.^2 - 8*x.*y + 2*x).*(2*z - 4*y - 4*x + 2*x.^2.*y.^2 + 8*x.*y - 4*x.*y.^2 - 4*x.^2.*y - 4*x.^2.*z + 2*x.^2 + 2*y.^2 + 16*x.^2.*y.*z - 12*x.^2.*y.^2.*z + 2))./w.^2).*(y>0.5);
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:,:) = (((8*y.^2.*z.^3)/5 + 2*x.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 8*y.^2)./w.^2 - (2*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2).*(x.*((8*y.^2.*z.^3)/5 + 8*y.^2) - (4*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8*y.^4 + 4*y.^2) + (8*y.^4.*z.^3)/5 - 4*y.^2) + 2))./w.^3).*(y<0.5) + ...
%!     ((32*y + 2*x.*(12*y.^2 - 16*y + 4) + z.^3.*((32*y)/5 + 2*x.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (24*y.^2)/5 - 8/5) - 24*y.^2 + 2*x.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) - 8)./w.^2 - (2*(z.^3.*(x.^2.*((72*y.^4)/5 - (192*y.^3)/5 + (176*y.^2)/5 - (64*y)/5 + 8/5) - (16*y)/5 - x.*((24*y.^2)/5 - (32*y)/5 + 8/5) + (12*y.^2)/5 + 4/5) - x.*(24*y.^2 - 32*y + 8) + x.^2.*(12*y.^2 - 16*y + 4) + x.^2.*z.^2.*(72*y.^4 - 192*y.^3 + 164*y.^2 - 48*y + 4) + 2).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3).*(y>0.5);
%! d2Fduu(2,:,:,:) = ((z.^3.*((24*x.*y.^4)/5 + (12*y.^2)/5) + (2*y.^2.*(20*y.^2 - 40*y + 20))/5 + z.^2.*((2*y.^2.*(40*y - 20*y.^2))/5 + 16*x.*y.^4))./w.^2 - (2*(z.^2.*(8*x.^2.*y.^4 - y.^2.*(8*y - 4*y.^2) + (2*x.*y.^2.*(40*y - 20*y.^2))/5) + z.^3.*((12*x.^2.*y.^4)/5 + (12*x.*y.^2)/5 - (6*y.^2)/5) + (2*x.*y.^2.*(20*y.^2 - 40*y + 20))/5).*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2))./w.^3).*(y<0.5) + ...
%!     (((2*(3*y.^2 - 4*y + 1).*(36*x.*y.^2 - 48*x.*y + 12*x - 6).*z.^3)/5 - (2*(3*y.^2 - 4*y + 1).*(160*x.*y - 40*x - 120*x.*y.^2 + 20*y.^2 + 10).*z.^2)/5 + (2*(20*y.^2 - 10).*(3*y.^2 - 4*y + 1))/5)./w.^2 - (2*((2*(3*y.^2 - 4*y + 1).*(18*x.^2.*y.^2 - 24*x.^2.*y + 6*x.^2 - 6*x + 3).*z.^3)/5 + (2*(3*y.^2 - 4*y + 1).*(60*x.^2.*y.^2 - 80*x.^2.*y + 20*x.^2 - 20*x.*y.^2 - 10*x + 10*y.^2 + 5).*z.^2)/5 - (2*(10*x - 20*x.*y.^2).*(3*y.^2 - 4*y + 1))/5).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3).*(y>0.5);
%! d2Fduu(3,:,:,:) =  (((16*x.*y.^4 + 8*y.^2).*z.^3 + (4*y.*(4*y.^2 - 6*y.^3) - 8*x.*y.*(2*y.^2 - 3*y.^3)).*z.^2 - 4*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2) + 8*x.*y.*(4*y.^2 - 6*y.^3))./w.^2 - (2*(2*y.^2.*z.^2 + 4*x.*y.^2 - 4*x.*y.^2.*z.^2).*(4*y.*(3*y - 2) + z.^3.*(8*x.^2.*y.^4 + 8*x.*y.^2 - 4*y.^2) - z.^2.*(4*y.*(2*y.^2 - 3*y.^3).*x.^2 - 4*y.*(4*y.^2 - 6*y.^3).*x + 4*y.*(2*y.^2 - 3*y.^3)) + 4*x.^2.*y.*(4*y.^2 - 6*y.^3) - 4*x.*y.*(- 6*y.^3 + 4*y.^2 + 3*y - 2)))./w.^3).*(y<0.5) + ...
%!  (-((4*(6*y - 2).*(y - 1) - 8*x.*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2)).*z.^3 + (4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2) - 8*x.*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)).*z.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) + 8*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2))./w.^2 - (2*(z.^2.*(4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1).*x.^2 - 4*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2).*x + 4*(y - 1).*(3*y.^3 - 7*y.^2 + 5*y - 1)) - 4*(y - 1).^2 + z.^3.*(4*(y - 1).*(18*y.^3 - 30*y.^2 + 14*y - 2).*x.^2 - 4*(6*y - 2).*(y - 1).*x + 4*(3*y - 1).*(y - 1)) + 4*x.*(y - 1).*(6*y.^3 - 14*y.^2 + 11*y - 3) - 4*x.^2.*(y - 1).*(6*y.^3 - 14*y.^2 + 10*y - 2)).*(4*x + 6*y.^2.*z.^2 - 16*x.*y + 12*x.*y.^2 - 4*x.*z.^2 - 8*y.*z.^2 + 2*z.^2 + 16*x.*y.*z.^2 - 12*x.*y.^2.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:,:) = ((((8.*x.^2.*(6.*z.^3 - 6.*z.^5))/5 + (8.*x.^4.*(10.*z.^4 - 20.*z.^2 + 10))/5 - (8.*x.^3.*(- 4.*z.^5 + 10.*z.^4 + 4.*z.^3 - 30.*z.^2 + 20))/5 + (16.*x.*z.^5)/5).*y.^3 + ((8.*x.*(2.*z.^3 - 10.*z.^2 + 10))/5 + (8.*x.^2.*(15.*z.^2 - 15))/5 - (8.*z.^3)/5).*y)./w.^3) .* (y<0.5) + ...
%!     (-(x.^4.*((8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10).*z.^4)/5 - (8.*(3.*y - 2).*(60.*y.^2 - 80.*y + 20).*z.^2)/5 + (8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10))/5) - x.^3.*(- (8.*(3.*y - 2).*(12.*y.^2 - 16.*y + 4).*z.^5)/5 + (8.*(3.*y - 2).*(30.*y.^2 - 40.*y + 10).*z.^4)/5 + (8.*(3.*y - 2).*(12.*y.^2 - 16.*y + 4).*z.^3)/5 - (8.*(3.*y - 2).*(90.*y.^2 - 120.*y + 30).*z.^2)/5 + (8.*(3.*y - 2).*(60.*y.^2 - 80.*y + 20))/5) + z.^3.*((24.*y)/5 - 16/5) - x.^2.*((8.*(3.*y - 2).*(18.*y.^2 - 24.*y + 6).*z.^5)/5 - (8.*(3.*y - 2).*(18.*y.^2 - 24.*y + 6).*z.^3)/5 + (72.*y - 48).*z.^2 - 72.*y + 48) + x.*((8.*(3.*y - 2).*(6.*y.^2 - 8.*y + 2).*z.^5)/5 + (32/5 - (48.*y)/5).*z.^3 + (48.*y - 32).*z.^2 - 48.*y + 32))./(-w).^3) .* (y>0.5);
%! d2Fduv(2,:,:,:) = ((((4.*x.^2.*(60.*z.^2 - 60.*z.^4))/5 + (4.*x.^3.*(40.*z.^4 - 80.*z.^2 + 40))/5 + 16.*x.*z.^4).*y.^4 + ((4.*x.^2.*(18.*z.^3 - 18.*z.^5))/5 + (4.*x.^3.*(12.*z.^5 - 12.*z.^3 + 40.*z.^2 - 40))/5 + (4.*x.*(6.*z.^5 - 40.*z.^2 + 40))/5 + 16.*z.^2).*y.^3 + ((4.*x.*(60.*z.^2 - 60))/5 - 24.*z.^2).*y.^2 + ((4.*x.*(6.*z.^3 + 20))/5 - (12.*z.^3)/5).*y)./w.^3) .* (y<0.5) + ...
%!     ((z.^3.*(((432.*y.^3)/5 - (864.*y.^2)/5 + (528.*y)/5 - 96/5).*x.^3 + (- (648.*y.^3)/5 + (1296.*y.^2)/5 - (792.*y)/5 + 144/5).*x.^2 + ((72.*y)/5 - 48/5).*x - (36.*y)/5 + 24/5) - x.^3.*(192.*y.^4 - 496.*y.^3 + 480.*y.^2 - 208.*y + 32) + z.^4.*((- 192.*y.^4 + 208.*y.^3 + 96.*y.^2 - 144.*y + 32).*x.^3 + (288.*y.^4 - 312.*y.^3 - 144.*y.^2 + 216.*y - 48).*x.^2 + (- 96.*y.^4 + 104.*y.^3 + 48.*y.^2 - 72.*y + 16).*x) + x.*(- 96.*y.^3 + 96.*y.^2 + 8.*y - 16) + z.^2.*(x.^2.*(- 288.*y.^4 + 312.*y.^3 + 144.*y.^2 - 216.*y + 48) - 20.*y - x.^3.*(- 384.*y.^4 + 704.*y.^3 - 384.*y.^2 + 64.*y) + x.*(96.*y.^3 - 96.*y.^2 + 40.*y - 16) + 48.*y.^2 - 48.*y.^3 + 8) - z.^5.*(((432.*y.^3)/5 - (864.*y.^2)/5 + (528.*y)/5 - 96/5).*x.^3 + (- (648.*y.^3)/5 + (1296.*y.^2)/5 - (792.*y)/5 + 144/5).*x.^2 + ((216.*y.^3)/5 - (432.*y.^2)/5 + (264.*y)/5 - 48/5).*x))./(-w).^3) .* (y>0.5);
%! d2Fduv(3,:,:,:) = (((x.^2.*(48.*z.^2 - 48.*z.^4) - x.^4.*(16.*z.^4 - 48.*z.^2 + 32) + x.^3.*(48.*z.^4 - 96.*z.^2 + 32) + 16.*x.*z.^4).*y.^4 + (x.^2.*(- 48.*z.^5 + 48.*z.^3 + 144.*z.^2 - 144) - x.^3.*(- 32.*z.^5 + 32.*z.^3 + 48.*z.^2 - 48) + x.*(16.*z.^5 - 144.*z.^2 + 96) + 48.*z.^2).*y.^3 + (x.*(96.*z.^2 - 48) + x.^3.*(48.*z.^2 - 48) - x.^2.*(120.*z.^2 - 96) - 24.*z.^2).*y.^2 + (x.*(16.*z.^3 - 24) - 8.*z.^3 + 24).*y + 8.*x - 8)./w.^3) .* (y<0.5) + ...
%!     ((8.*y - x.^4.*(96.*y.^4 - 320.*y.^3 + 384.*y.^2 - 192.*y + 32) + x.^3.*(96.*y.^4 - 368.*y.^3 + 528.*y.^2 - 336.*y + 80) + z.^3.*((288.*y.^3 - 576.*y.^2 + 352.*y - 64).*x.^3 + (- 432.*y.^3 + 864.*y.^2 - 528.*y + 96).*x.^2 + (48.*y - 32).*x - 24.*y + 16) - x.*(96.*y.^3 - 240.*y.^2 + 200.*y - 56) - z.^4.*((48.*y.^4 - 160.*y.^3 + 192.*y.^2 - 96.*y + 16).*x.^4 + (- 144.*y.^4 + 480.*y.^3 - 576.*y.^2 + 288.*y - 48).*x.^3 + (144.*y.^4 - 480.*y.^3 + 576.*y.^2 - 288.*y + 48).*x.^2 + (- 48.*y.^4 + 160.*y.^3 - 192.*y.^2 + 96.*y - 16).*x) + z.^2.*(x.^4.*(144.*y.^4 - 480.*y.^3 + 576.*y.^2 - 288.*y + 48) - 96.*y + x.^2.*(144.*y.^4 - 624.*y.^3 + 984.*y.^2 - 672.*y + 168) - x.^3.*(288.*y.^4 - 1008.*y.^3 + 1296.*y.^2 - 720.*y + 144) + x.*(144.*y.^3 - 384.*y.^2 + 336.*y - 96) + 120.*y.^2 - 48.*y.^3 + 24) - z.^5.*((288.*y.^3 - 576.*y.^2 + 352.*y - 64).*x.^3 + (- 432.*y.^3 + 864.*y.^2 - 528.*y + 96).*x.^2 + (144.*y.^3 - 288.*y.^2 + 176.*y - 32).*x) + x.^2.*(144.*y.^3 - 384.*y.^2 + 336.*y - 96) - 8)./(-w).^3) .* (y>0.5);
%! d2Fduw = zeros ([3, size(x)]);
%! d2Fduw(1,:,:,:) = ((x.^2.*((24.*y.^4.*z.^2)/5 + 2.*z.*(8.*y.^4 + 4.*y.^2)) - (12.*y.^2.*z.^2)/5 + (24.*x.*y.^2.*z.^2)/5)./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(x.*((8.*y.^2.*z.^3)/5 + 8.*y.^2) - (4.*y.^2.*z.^3)/5 + x.^2.*(z.^2.*(8.*y.^4 + 4.*y.^2) + (8.*y.^4.*z.^3)/5 - 4.*y.^2) + 2))./w.^3) .* (y<0.5) + ...
%!     (-((- (4.*(3.*y - 1).*(y - 1).*(36.*y.^4 - 96.*y.^3 + 88.*y.^2 - 32.*y + 4).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(36.*y.^4 - 96.*y.^3 + 100.*y.^2 - 48.*y + 8).*x.^3)/5 - (4.*(3.*y - 1).*(y - 1).*(18.*y.^2 - 24.*y + 6).*x.^2)/5 + (4.*(3.*y - 1).*(y - 1).*(6.*y.^2 - 8.*y + 2).*x)/5).*z.^4 + ((4.*x.^3.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 820.*y.^2 - 240.*y + 20))/5 - (4.*x.^4.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 820.*y.^2 - 240.*y + 20))/5).*z.^3 + (- (4.*(3.*y - 1).*(y - 1).*(108.*y.^4 - 288.*y.^3 + 264.*y.^2 - 96.*y + 12).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(36.*y.^2 - 48.*y + 12).*x.^3)/5 - (24.*(3.*y - 1).*(y - 1).*x)/5 + (12.*(3.*y - 1).*(y - 1))/5).*z.^2 + (- (4.*(3.*y - 1).*(y - 1).*(360.*y.^4 - 960.*y.^3 + 940.*y.^2 - 400.*y + 60).*x.^4)/5 + (4.*(3.*y - 1).*(y - 1).*(360.*y.^2 - 480.*y + 120).*x.^3)/5 - (4.*(3.*y - 1).*(y - 1).*(180.*y.^2 - 240.*y + 90).*x.^2)/5 + 16.*(3.*y - 1).*(y - 1).*x).*z)./(-w).^3) .* (y>0.5);
%! d2Fduw(2,:,:,:) = ((2.*z.*(8.*x.^2.*y.^4 - y.^2.*(8.*y - 4.*y.^2) + (2.*x.*y.^2.*(40.*y - 20.*y.^2))/5) + 3.*z.^2.*((12.*x.^2.*y.^4)/5 + (12.*x.*y.^2)/5 - (6.*y.^2)/5))./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(z.^2.*(8.*x.^2.*y.^4 - y.^2.*(8.*y - 4.*y.^2) + (2.*x.*y.^2.*(40.*y - 20.*y.^2))/5) + z.^3.*((12.*x.^2.*y.^4)/5 + (12.*x.*y.^2)/5 - (6.*y.^2)/5) + (2.*x.*y.^2.*(20.*y.^2 - 40.*y + 20))/5))./w.^3) .* (y<0.5) + ...
%!     (((6.*(3.*y.^2 - 4.*y + 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 6.*x + 3).*z.^2)/5 + (4.*(3.*y.^2 - 4.*y + 1).*(60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2 - 20.*x.*y.^2 - 10.*x + 10.*y.^2 + 5).*z)/5)./w.^2 - (2.*((2.*(3.*y.^2 - 4.*y + 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 6.*x + 3).*z.^3)/5 + (2.*(3.*y.^2 - 4.*y + 1).*(60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2 - 20.*x.*y.^2 - 10.*x + 10.*y.^2 + 5).*z.^2)/5 - (2.*(10.*x - 20.*x.*y.^2).*(3.*y.^2 - 4.*y + 1))/5).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fduw(3,:,:,:) = (- (2.*z.*(4.*y.*(2.*y.^2 - 3.*y.^3).*x.^2 - 4.*y.*(4.*y.^2 - 6.*y.^3).*x + 4.*y.*(2.*y.^2 - 3.*y.^3)) - 3.*z.^2.*(8.*x.^2.*y.^4 + 8.*x.*y.^2 - 4.*y.^2))./w.^2 - (2.*(4.*x.*y.^2.*z - 4.*x.^2.*y.^2.*z).*(4.*y.*(3.*y - 2) + z.^3.*(8.*x.^2.*y.^4 + 8.*x.*y.^2 - 4.*y.^2) - z.^2.*(4.*y.*(2.*y.^2 - 3.*y.^3).*x.^2 - 4.*y.*(4.*y.^2 - 6.*y.^3).*x + 4.*y.*(2.*y.^2 - 3.*y.^3)) + 4.*x.^2.*y.*(4.*y.^2 - 6.*y.^3) - 4.*x.*y.*(- 6.*y.^3 + 4.*y.^2 + 3.*y - 2)))./w.^3) .* (y<0.5) + ...
%!     ((2.*z.*(4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1).*x.^2 - 4.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2).*x + 4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1)) + 3.*z.^2.*(4.*(y - 1).*(18.*y.^3 - 30.*y.^2 + 14.*y - 2).*x.^2 - 4.*(6.*y - 2).*(y - 1).*x + 4.*(3.*y - 1).*(y - 1)))./w.^2 - (2.*(z.^2.*(4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1).*x.^2 - 4.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2).*x + 4.*(y - 1).*(3.*y.^3 - 7.*y.^2 + 5.*y - 1)) - 4.*(y - 1).^2 + z.^3.*(4.*(y - 1).*(18.*y.^3 - 30.*y.^2 + 14.*y - 2).*x.^2 - 4.*(6.*y - 2).*(y - 1).*x + 4.*(3.*y - 1).*(y - 1)) + 4.*x.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 11.*y - 3) - 4.*x.^2.*(y - 1).*(6.*y.^3 - 14.*y.^2 + 10.*y - 2)).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:,:) = (-(8.*x.*(x - 1).*(z.^3 + 5.*x.*z.^2 - 5.*x).*(- 6.*x.^2.*y.^2.*z.^2 + 6.*x.^2.*y.^2 + 6.*x.*y.^2.*z.^2 - 1))/5./w.^3) .* (y<0.5) + ...
%!     ((8.*x.*(x - 1).*(z.^3 + 5.*x.*z.^2 - 5.*x).*(- 54.*x.^2.*y.^2.*z.^2 + 54.*x.^2.*y.^2 + 72.*x.^2.*y.*z.^2 - 72.*x.^2.*y - 26.*x.^2.*z.^2 + 26.*x.^2 + 54.*x.*y.^2.*z.^2 - 72.*x.*y.*z.^2 + 26.*x.*z.^2 + 3))/5./(-w).^3) .* (y>0.5);
%! d2Fdvv(2,:,:,:) = ((2.*((8.*x.*z.^2 - x.^2.*(8.*z.^2 - 8)).*y.^2 + ((12.*x.*z.^3)/5 - x.^2.*((12.*z.^3)/5 + 8) + 4).*y - 4).*(- 4.*y.*x.^2.*z.^2 + 4.*y.*x.^2 + 4.*y.*x.*z.^2))./w.^3 - ((12.*x.*z.^3)/5 + 2.*y.*(8.*x.*z.^2 - x.^2.*(8.*z.^2 - 8)) - x.^2.*((12.*z.^3)/5 + 8) + 4)./w.^2) .* (y<0.5) + ...
%!     ((z.^2.*(x.*(32.*y + 4) - x.^2.*(32.*y + 4)) + x.^2.*(32.*y - 20) + z.^3.*((36.*x)/5 - (36.*x.^2)/5) + 4)./w.^2 - (2.*(4.*y + z.^3.*(x.*((36.*y)/5 - 24/5) - x.^2.*((36.*y)/5 - 24/5)) + z.^2.*(x.*(16.*y.^2 + 4.*y - 8) - x.^2.*(16.*y.^2 + 4.*y - 8)) + x.^2.*(16.*y.^2 - 20.*y + 8)).*(8.*x.^2.*z.^2 + 12.*x.^2.*y - 8.*x.*z.^2 - 8.*x.^2 + 12.*x.*y.*z.^2 - 12.*x.^2.*y.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fdvv(3,:,:,:) = ((2.*y.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(3.*x - 3).*(x - 1) + 8.*x.*z.^3.*(x - 1))./w.^2 - (2.*(4.*(x - 1).^2 - y.*(4.*(3.*x - 3).*(x - 1) - 8.*x.*z.^3.*(x - 1)) + y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1))).*(- 4.*y.*x.^2.*z.^2 + 4.*y.*x.^2 + 4.*y.*x.*z.^2))./w.^3) .* (y<0.5) + ...
%!     ((4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1) + 2.*y.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 24.*x.*z.^3.*(x - 1) - 4.*z.^2.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x))./w.^2 - (2.*(y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(x - 1).*(2.*x.^3 - 2.*x.^2 + x - 1) - y.*(24.*x.*(x - 1).*z.^3 + 4.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z.^2 - 4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1)) + 16.*x.*z.^3.*(x - 1) + 4.*z.^2.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x)).*(8.*x.^2.*z.^2 + 12.*x.^2.*y - 8.*x.*z.^2 - 8.*x.^2 + 12.*x.*y.*z.^2 - 12.*x.^2.*y.*z.^2))./(-w).^3) .* (y>0.5);
%! d2Fdvw = zeros ([3, size(x)]);
%! d2Fdvw(1,:,:,:) = (((8.*x.*z.*(x - 1).*(20.*x.^3.*z.^2 - 20.*x.^3 + 2.*x.^2.*z.^3 - 20.*x.^2.*z.^2 + 6.*x.^2.*z + 40.*x.^2 - 2.*x.*z.^3).*y.^3)/5 + (8.*x.*z.*(10.*x + 3.*z).*(x - 1).*y)/5)./w.^3) .* (y<0.5) + ...
%!     (((8.*x.*(3.*y - 2).*(x - 1).*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*z.^4)/5 + (8.*x.*(3.*y - 2).*(x - 1).*(- 60.*x.^3.*y.^2 + 80.*x.^3.*y - 20.*x.^3 + 60.*x.^2.*y.^2 - 80.*x.^2.*y + 20.*x.^2).*z.^3)/5 - (8.*x.*(3.*y - 2).*(x - 1).*(18.*x.^2.*y.^2 - 24.*x.^2.*y + 6.*x.^2 - 3).*z.^2)/5 + (8.*x.*(3.*y - 2).*(x - 1).*(60.*x.^3.*y.^2 - 80.*x.^3.*y + 20.*x.^3 - 120.*x.^2.*y.^2 + 160.*x.^2.*y - 40.*x.^2 + 10.*x).*z)/5)./(-w).^3) .* (y>0.5);
%! d2Fdvw(2,:,:,:) = ((4.*x.*y.*z.*(x - 1).*(40.*x.^2.*y.^3.*z.^2 - 40.*x.^2.*y.^3 + 6.*x.^2.*y.^2.*z.^3 + 18.*x.^2.*y.^2.*z + 80.*x.^2.*y.^2 - 40.*x.*y.^3.*z.^2 - 6.*x.*y.^2.*z.^3 - 40.*y.^2 + 60.*y + 9.*z))/5./w.^3) .* (y<0.5) + ...
%!     (-((4.*x.*(x - 1).*(54.*x.^2.*y.^3 - 108.*x.^2.*y.^2 + 66.*x.^2.*y - 12.*x.^2 - 54.*x.*y.^3 + 108.*x.*y.^2 - 66.*x.*y + 12.*x).*z.^4)/5 + (4.*x.*(x - 1).*(240.*x.^2.*y.^4 - 260.*x.^2.*y.^3 - 120.*x.^2.*y.^2 + 180.*x.^2.*y - 40.*x.^2 - 240.*x.*y.^4 + 260.*x.*y.^3 + 120.*x.*y.^2 - 180.*x.*y + 40.*x).*z.^3)/5 - (4.*x.*(x - 1).*(- 162.*x.^2.*y.^3 + 324.*x.^2.*y.^2 - 198.*x.^2.*y + 36.*x.^2 + 27.*y - 18).*z.^2)/5 - (4.*x.*(x - 1).*(240.*x.^2.*y.^4 - 980.*x.^2.*y.^3 + 1320.*x.^2.*y.^2 - 700.*x.^2.*y + 120.*x.^2 + 120.*y.^3 - 120.*y.^2 + 50.*y - 20).*z)/5)./(-w).^3) .* (y>0.5);
%! d2Fdvw(3,:,:,:) = (-(y.^3.*(8.*x.*z.*(x - 1).*(12.*x.^2 - 24.*x + 12) - 48.*x.^3.*z.^2.*(x - 1) + 8.*x.*z.^4.*(2.*x - 2.*x.^2).*(x - 1)) + y.^4.*(8.*x.*(x - 1).*(- 4.*x.^4 + 12.*x.^3 - 12.*x.^2 + 4.*x).*z.^3 + 8.*x.*(x - 1).*(4.*x.^4 - 8.*x.^3 + 4.*x.^2).*z) - 24.*x.*y.*z.^2.*(x - 1) - 8.*x.*y.^2.*z.*(x - 1).*(6.*x.^2 - 12.*x + 6))./w.^3) .* (y<0.5) + ...
%!     ((8.*z.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x) - y.*(72.*x.*(x - 1).*z.^2 + 8.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z) + 48.*x.*z.^2.*(x - 1) + 8.*y.^2.*z.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x))./w.^2 - (2.*(y.^2.*(4.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x).*z.^2 + 4.*(2.*x.^2 - 2.*x.^3).*(x - 1)) - 4.*(x - 1).*(2.*x.^3 - 2.*x.^2 + x - 1) - y.*(24.*x.*(x - 1).*z.^3 + 4.*(x - 1).*(4.*x.^3 - 8.*x.^2 + 4.*x).*z.^2 - 4.*(x - 1).*(4.*x.^3 - 4.*x.^2 + x - 1)) + 16.*x.*z.^3.*(x - 1) + 4.*z.^2.*(x - 1).*(2.*x.^3 - 4.*x.^2 + 2.*x)).*(- 12.*z.*x.^2.*y.^2 + 16.*z.*x.^2.*y - 4.*z.*x.^2 + 12.*z.*x.*y.^2 - 16.*z.*x.*y + 4.*z.*x))./(-w).^3) .* (y>0.5);
%! d2Fdww = zeros ([3, size(x)]);
%! d2Fdww(1,:,:,:) = ((32.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(5.*x + z + 10.*x.^2.*y.^2 + 2.*x.^2.*y.^2.*z))./(5.*w.^3) - (8.*x.*y.^2.*(x - 1).*(15.*x + z + 30.*x.^2.*y.^2 + 2.*x.^2.*y.^2.*z))/5./w.^2) .* (y<0.5) + ...
%!     (((8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(36.*x.^4.*y.^4 - 96.*x.^4.*y.^3 + 88.*x.^4.*y.^2 - 32.*x.^4.*y + 4.*x.^4 - 36.*x.^3.*y.^4 + 96.*x.^3.*y.^3 - 88.*x.^3.*y.^2 + 32.*x.^3.*y - 4.*x.^3 - 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*z.^3)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(540.*x.^4.*y.^4 - 1440.*x.^4.*y.^3 + 1320.*x.^4.*y.^2 - 480.*x.^4.*y + 60.*x.^4 - 540.*x.^3.*y.^4 + 1440.*x.^3.*y.^3 - 1410.*x.^3.*y.^2 + 600.*x.^3.*y - 90.*x.^3 + 90.*x.^2.*y.^2 - 120.*x.^2.*y + 30.*x.^2).*z.^2)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(108.*x.^4.*y.^4 - 288.*x.^4.*y.^3 + 264.*x.^4.*y.^2 - 96.*x.^4.*y + 12.*x.^4 - 36.*x.^2.*y.^2 + 48.*x.^2.*y - 12.*x.^2 + 3).*z)/5 + (8.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(180.*x.^4.*y.^4 - 480.*x.^4.*y.^3 + 440.*x.^4.*y.^2 - 160.*x.^4.*y + 20.*x.^4 - 30.*x.^3.*y.^2 + 40.*x.^3.*y - 10.*x.^3 - 30.*x.^2.*y.^2 + 40.*x.^2.*y - 10.*x.^2 + 5.*x))/5)./(-w).^3) .* (y>0.5);
%! d2Fdww(2,:,:,:) = ((16.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(20.*y + 3.*z + 20.*x.^2.*y.^2 - 10.*y.^2 + 6.*x.^2.*y.^2.*z))./(5.*w.^3) - (12.*x.*y.^2.*(x - 1).*(20.*y + z + 20.*x.^2.*y.^2 - 10.*y.^2 + 2.*x.^2.*y.^2.*z))/5./w.^2) .* (y<0.5) + ...
%!     (((4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(108.*x.^4.*y.^4 - 288.*x.^4.*y.^3 + 264.*x.^4.*y.^2 - 96.*x.^4.*y + 12.*x.^4 - 108.*x.^3.*y.^4 + 288.*x.^3.*y.^3 - 264.*x.^3.*y.^2 + 96.*x.^3.*y - 12.*x.^3 - 18.*x.^2.*y.^2 + 24.*x.^2.*y - 6.*x.^2 + 18.*x.*y.^2 - 24.*x.*y + 6.*x).*z.^3)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(1080.*x.^4.*y.^4 - 2880.*x.^4.*y.^3 + 2640.*x.^4.*y.^2 - 960.*x.^4.*y + 120.*x.^4 - 1080.*x.^3.*y.^4 + 2880.*x.^3.*y.^3 - 2640.*x.^3.*y.^2 + 960.*x.^3.*y - 120.*x.^3 - 180.*x.^2.*y.^4 + 240.*x.^2.*y.^3 - 150.*x.^2.*y.^2 + 120.*x.^2.*y - 30.*x.^2 + 180.*x.*y.^4 - 240.*x.*y.^3 + 150.*x.*y.^2 - 120.*x.*y + 30.*x).*z.^2)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(324.*x.^4.*y.^4 - 864.*x.^4.*y.^3 + 792.*x.^4.*y.^2 - 288.*x.^4.*y + 36.*x.^4 - 108.*x.^2.*y.^2 + 144.*x.^2.*y - 36.*x.^2 + 9).*z)/5 + (4.*x.*(3.*y - 1).*(x - 1).*(y - 1).*(360.*x.^4.*y.^4 - 960.*x.^4.*y.^3 + 880.*x.^4.*y.^2 - 320.*x.^4.*y + 40.*x.^4 - 60.*x.^2.*y.^4 + 80.*x.^2.*y.^3 - 110.*x.^2.*y.^2 + 120.*x.^2.*y - 30.*x.^2 + 10.*y.^2 + 5))/5)./(-w).^3) .* (y>0.5);
%! d2Fdww(3,:,:,:) = ((32.*x.*y.^2.*(2.*x.^2.*y.^2 + 1).*(x - 1).*(2.*y + z - 3.*x.^2.*y.^2 - 4.*x.*y + 6.*x.*y.^2 + 2.*x.^2.*y - 3.*y.^2 + 2.*x.^2.*y.^2.*z))./w.^3 - (8.*x.*y.^2.*(x - 1).*(6.*y + z - 9.*x.^2.*y.^2 - 12.*x.*y + 18.*x.*y.^2 + 6.*x.^2.*y - 9.*y.^2 + 2.*x.^2.*y.^2.*z))./w.^2) .* (y<0.5) + ...
%!    ((2.*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).*(2.*z - 4.*y - 4.*x + 2.*x.^2.*y.^2 + 8.*x.*y - 4.*x.*y.^2 - 4.*x.^2.*y - 4.*x.^2.*z + 2.*x.^2 + 2.*y.^2 + 16.*x.^2.*y.*z - 12.*x.^2.*y.^2.*z + 2))./w.^2 - (8.*z.^2.*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x).^2.*(2.*z - 4.*y - 4.*x + 2.*x.^2.*y.^2 + 8.*x.*y - 4.*x.*y.^2 - 4.*x.^2.*y - 4.*x.^2.*z + 2.*x.^2 + 2.*y.^2 + 16.*x.^2.*y.*z - 12.*x.^2.*y.^2.*z + 2))./(-w).^3 - (4.*z.*(12.*x.^2.*y.^2 - 16.*x.^2.*y + 4.*x.^2 - 2).*(- 6.*x.^2.*y.^2 + 8.*x.^2.*y - 2.*x.^2 + 6.*x.*y.^2 - 8.*x.*y + 2.*x))./w.^2) .* (y>0.5);
%! assert (F, pnt, 1e3*eps)
%! assert (dFdu, jac{1}, 1e3*eps)
%! assert (dFdv, jac{2}, 1e3*eps)
%! assert (dFdw, jac{3}, 1e3*eps)
%! assert (d2Fduu, hess{1,1}, 1e3*eps)
%! assert (d2Fduv, hess{1,2}, 1e3*eps)
%! assert (d2Fduw, hess{1,3}, 1e3*eps)
%! assert (d2Fduv, hess{2,1}, 1e3*eps)
%! assert (d2Fdvv, hess{2,2}, 1e3*eps)
%! assert (d2Fdvw, hess{2,3}, 1e3*eps)
%! assert (d2Fduw, hess{3,1}, 1e3*eps)
%! assert (d2Fdvw, hess{3,2}, 1e3*eps)
%! assert (d2Fdww, hess{3,3}, 1e3*eps)



%!test
%! nrb = nrbextrude (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [0 0 1]);
%! nrb = nrbdegelev (nrb, [1 1 1]);
%! nrb.coefs (4,2,2,2) = 1.1;
%! [dnrb, dnrb2] = nrbderiv (nrb);
%! X = linspace (0, 1, 24); Y = linspace (0, 1, 24); Z = linspace (0, 1, 24);
%! [pnt, jac, hess] = nrbdeval (nrb, dnrb, dnrb2, {X Y Z});
%! [y, x, z] = meshgrid (X, Y, Z);
%! F = zeros ([3, size(x)]);
%! F(1,:,:,:) = (5.*x)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! F(2,:,:,:) = (5.*y)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! F(3,:,:,:) = (5.*z)./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5);
%! dFdu = zeros ([3, size(x)]);
%! dFdu(1,:,:,:) = ((z.*(20.*y - 20.*y.^2) - z.^2.*(20.*y - 20.*y.^2)).*x.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! dFdu(2,:,:,:) = (y.^2.*(5.*z.*(8.*x - 4) - 5.*z.^2.*(8.*x - 4)) - y.^3.*(5.*z.*(8.*x - 4) - 5.*z.^2.*(8.*x - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdu(3,:,:,:) = (z.^2.*(5.*y.*(8.*x - 4) - 5.*y.^2.*(8.*x - 4)) - z.^3.*(5.*y.*(8.*x - 4) - 5.*y.^2.*(8.*x - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdv = zeros ([3, size(x)]);
%! dFdv(1,:,:,:) = (x.^2.*(5.*z.*(8.*y - 4) - 5.*z.^2.*(8.*y - 4)) - x.^3.*(5.*z.*(8.*y - 4) - 5.*z.^2.*(8.*y - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdv(2,:,:,:) = ((z.*(20.*x - 20.*x.^2) - z.^2.*(20.*x - 20.*x.^2)).*y.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! dFdv(3,:,:,:) = (z.^2.*(5.*x.*(8.*y - 4) - 5.*x.^2.*(8.*y - 4)) - z.^3.*(5.*x.*(8.*y - 4) - 5.*x.^2.*(8.*y - 4)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw = zeros ([3, size(x)]);
%! dFdw(1,:,:,:) = (x.^2.*(y.*(40.*z - 20) - y.^2.*(40.*z - 20)) - x.^3.*(y.*(40.*z - 20) - y.^2.*(40.*z - 20)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw(2,:,:,:) = (y.^2.*(x.*(40.*z - 20) - x.^2.*(40.*z - 20)) - y.^3.*(x.*(40.*z - 20) - x.^2.*(40.*z - 20)))./((- 4.*x.^2.*y.^2 + 4.*x.^2.*y + 4.*x.*y.^2 - 4.*x.*y).*z.^2 + (4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y).*z + 5).^2;
%! dFdw(3,:,:,:) = ((y.*(20.*x - 20.*x.^2) - y.^2.*(20.*x - 20.*x.^2)).*z.^2 + 25)./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^2;
%! d2Fduu = zeros ([3, size(x)]);
%! d2Fduu(1,:,:,:) = (40.*y.*z.*(y - 1).*(z - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z + 15.*x - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduu(2,:,:,:) = (40.*y.^2.*z.*(y - 1).*(z - 1).*(4.*y.^2.*z.^2 - 4.*y.^2.*z - 4.*y.*z.^2 + 4.*y.*z + 5) - 40.*x.*y.^2.*z.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z) + 40.*x.^2.*y.^2.*z.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduu(3,:,:,:) = (40.*y.*z.^2.*(y - 1).*(z - 1).*(4.*y.^2.*z.^2 - 4.*y.^2.*z - 4.*y.*z.^2 + 4.*y.*z + 5) - 40.*x.*y.*z.^2.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z) + 40.*x.^2.*y.*z.^2.*(y - 1).*(z - 1).*(12.*y.^2.*z.^2 - 12.*y.^2.*z - 12.*y.*z.^2 + 12.*y.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv = zeros ([3, size(x)]);
%! d2Fduv(1,:,:,:) = (20.*x.*z.*(2.*y - 1).*(z - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 15.*x - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv(2,:,:,:) = (20.*y.*z.*(2.*x - 1).*(z - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z + 15.*y - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduv(3,:,:,:) = (20.*z.^2.*(2.*x - 1).*(2.*y - 1).*(z - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw = zeros ([3, size(x)]);
%! d2Fduw(1,:,:,:) = (20.*x.*y.*(2.*z - 1).*(y - 1).*(4.*x.^3.*y.^2.*z.^2 - 4.*x.^3.*y.^2.*z - 4.*x.^3.*y.*z.^2 + 4.*x.^3.*y.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 15.*x - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw(2,:,:,:) = (20.*y.^2.*(2.*x - 1).*(2.*z - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fduw(3,:,:,:) = (20.*y.*z.*(2.*x - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.*z.^3 + 4.*x.^2.*y.*z.^2 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.^2.*z.^2 + 4.*x.*y.*z.^3 - 4.*x.*y.*z.^2 + 15.*z - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv = zeros ([3, size(x)]);
%! d2Fdvv(1,:,:,:) = (40.*x.^2.*z.*(x - 1).*(z - 1).*(4.*x.^2.*z.^2 - 4.*x.^2.*z - 4.*x.*z.^2 + 4.*x.*z + 5) + 40.*x.^2.*y.^2.*z.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z) - 40.*x.^2.*y.*z.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv(2,:,:,:) = (40.*x.*z.*(x - 1).*(z - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 15.*y - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvv(3,:,:,:) = (40.*x.*z.^2.*(x - 1).*(z - 1).*(4.*x.^2.*z.^2 - 4.*x.^2.*z - 4.*x.*z.^2 + 4.*x.*z + 5) + 40.*x.*y.^2.*z.^2.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z) - 40.*x.*y.*z.^2.*(x - 1).*(z - 1).*(12.*x.^2.*z.^2 - 12.*x.^2.*z - 12.*x.*z.^2 + 12.*x.*z))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw = zeros ([3, size(x)]);
%! d2Fdvw(1,:,:,:) = (20.*x.^2.*(2.*y - 1).*(2.*z - 1).*(x - 1).*(4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.^2.*z - 4.*x.^2.*y.*z.^2 + 4.*x.^2.*y.*z - 4.*x.*y.^2.*z.^2 + 4.*x.*y.^2.*z + 4.*x.*y.*z.^2 - 4.*x.*y.*z + 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw(2,:,:,:) = (20.*x.*y.*(2.*z - 1).*(x - 1).*(4.*x.^2.*y.^3.*z.^2 - 4.*x.^2.*y.^3.*z - 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z - 4.*x.*y.^3.*z.^2 + 4.*x.*y.^3.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z + 15.*y - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdvw(3,:,:,:) = (20.*x.*z.*(2.*y - 1).*(x - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.^2.*z.^2 - 4.*x.^2.*y.*z.^3 + 4.*x.^2.*y.*z.^2 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.^2.*z.^2 + 4.*x.*y.*z.^3 - 4.*x.*y.*z.^2 + 15.*z - 10))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww = zeros ([3, size(x)]);
%! d2Fdww(1,:,:,:) = (40.*x.^2.*y.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y + 5) + 40.*x.^2.*y.*z.^2.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y) - 40.*x.^2.*y.*z.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww(2,:,:,:) = (40.*x.*y.^2.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2 - 4.*x.^2.*y - 4.*x.*y.^2 + 4.*x.*y + 5) + 40.*x.*y.^2.*z.^2.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y) - 40.*x.*y.^2.*z.*(x - 1).*(y - 1).*(12.*x.^2.*y.^2 - 12.*x.^2.*y - 12.*x.*y.^2 + 12.*x.*y))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;
%! d2Fdww(3,:,:,:) = (40.*x.*y.*(x - 1).*(y - 1).*(4.*x.^2.*y.^2.*z.^3 - 4.*x.^2.*y.*z.^3 - 4.*x.*y.^2.*z.^3 + 4.*x.*y.*z.^3 + 15.*z - 5))./(- 4.*x.^2.*y.^2.*z.^2 + 4.*x.^2.*y.^2.*z + 4.*x.^2.*y.*z.^2 - 4.*x.^2.*y.*z + 4.*x.*y.^2.*z.^2 - 4.*x.*y.^2.*z - 4.*x.*y.*z.^2 + 4.*x.*y.*z + 5).^3;

