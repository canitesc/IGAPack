function varargout = nrbdeval (nurbs, dnurbs, varargin)

% NRBDEVAL: Evaluation of the derivative and second derivatives of NURBS curve, surface or volume.
%
%     [pnt, jac] = nrbdeval (crv, dcrv, tt)
%     [pnt, jac] = nrbdeval (srf, dsrf, {tu tv})
%     [pnt, jac] = nrbdeval (vol, dvol, {tu tv tw})
%     [pnt, jac, hess] = nrbdeval (crv, dcrv, dcrv2, tt)
%     [pnt, jac, hess] = nrbdeval (srf, dsrf, dsrf2, {tu tv})
%     [pnt, jac, hess] = nrbdeval (vol, dvol, dvol2, {tu tv tw})
%     [pnt, jac, hess, third_der] = nrbdeval (crv, dcrv, dcrv2, dcrv3, tt)
%     [pnt, jac, hess, third_der] = nrbdeval (srf, dsrf, dsrf2, dsrf3, {tu tv})
%     [pnt, jac, hess, third_der, fourth_der] = nrbdeval (crv, dcrv, dcrv2, dcrv3, dcrv4, tt)
%     [pnt, jac, hess, third_der, fourth_der] = nrbdeval (srf, dsrf, dsrf2, dsrf3, dsrf4, {tu tv})
%
% INPUTS:
%
%   crv,   srf,   vol   - original NURBS curve, surface or volume.
%   dcrv,  dsrf,  dvol  - NURBS derivative representation of crv, srf 
%                          or vol (see nrbderiv2)
%   dcrv2, dsrf2, dvol2 - NURBS second derivative representation of crv,
%                          srf or vol (see nrbderiv2)
%   dcrv3, dsrf3, dvol3 - NURBS third derivative representation of crv,
%                          srf or vol (see nrbderiv)
%   dcrv4, dsrf4, dvol4 - NURBS fourth derivative representation of crv,
%                          srf or vol (see nrbderiv)
%
%   tt     - parametric evaluation points
%            If the nurbs is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnt  - evaluated points.
%   jac  - evaluated first derivatives (Jacobian).
%   hess - evaluated second derivatives (Hessian).
%   third_der - evaluated third derivatives
%   fourth_der - evaluated fourth derivatives
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

if (nargin == 3)
  tt = varargin{1};
elseif (nargin == 4)
  dnurbs2 = varargin{1};
  tt = varargin{2};
elseif (nargin == 5)
  dnurbs2 = varargin{1};
  dnurbs3 = varargin{2};       
  tt = varargin{3}; 
elseif (nargin == 6)
  dnurbs2 = varargin{1};
  dnurbs3 = varargin{2};
  dnurbs4 = varargin{3};       
  tt = varargin{4}; 
else 
  error ('nrbrdeval: wrong number of input parameters')
end

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

[cp,cw] = nrbeval(nurbs, tt);

if (iscell(nurbs.knots))
  if (size(nurbs.knots,2) == 3)
  % NURBS structure represents a volume
    temp = cw(ones(3,1),:,:,:);
    pnt = cp./temp;
  
    [cup,cuw] = nrbeval (dnurbs{1}, tt);
    tempu = cuw(ones(3,1),:,:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp,cvw] = nrbeval (dnurbs{2}, tt);
    tempv = cvw(ones(3,1),:,:,:);
    jac{2} = (cvp-tempv.*pnt)./temp;

    [cwp,cww] = nrbeval (dnurbs{3}, tt);
    tempw = cww(ones(3,1),:,:,:);
    jac{3} = (cwp-tempw.*pnt)./temp;

% second derivatives
    if (nargout >= 3)
      if (exist ('dnurbs2'))
        [cuup, cuuw] = nrbeval (dnurbs2{1,1}, tt);
        tempuu = cuuw(ones(3,1),:,:,:);
        hess{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;
        clear cuup cuuw tempuu

        [cvvp, cvvw] = nrbeval (dnurbs2{2,2}, tt);
        tempvv = cvvw(ones(3,1),:,:,:);
        hess{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;
        clear cvvp cvvw tempvv

        [cwwp, cwww] = nrbeval (dnurbs2{3,3}, tt);
        tempww = cwww(ones(3,1),:,:,:);
        hess{3,3} = (cwwp - (2*cwp.*tempw + cp.*tempww)./temp + 2*cp.*tempw.^2./temp.^2)./temp;
        clear cwwp cwww tempww

        [cuvp, cuvw] = nrbeval (dnurbs2{1,2}, tt);
        tempuv = cuvw(ones(3,1),:,:,:);
        hess{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hess{2,1} = hess{1,2};
        clear cuvp cuvw tempuv

        [cuwp, cuww] = nrbeval (dnurbs2{1,3}, tt);
        tempuw = cuww(ones(3,1),:,:,:);
        hess{1,3} = (cuwp - (cup.*tempw + cwp.*tempu + cp.*tempuw)./temp + 2*cp.*tempu.*tempw./temp.^2)./temp;
        hess{3,1} = hess{1,3};
        clear cuwp cuww tempuw

        [cvwp, cvww] = nrbeval (dnurbs2{2,3}, tt);
        tempvw = cvww(ones(3,1),:,:,:);
        hess{2,3} = (cvwp - (cvp.*tempw + cwp.*tempv + cp.*tempvw)./temp + 2*cp.*tempv.*tempw./temp.^2)./temp;
        hess{3,2} = hess{2,3};
        clear cvwp cvww tempvw
      else
        warning ('nrbdeval: dnurbs2 missing. The second derivative is not computed');
        hess = [];
      end
      if (nargout > 3)
        warning ('nrbdeval: 3rd and 4th order derivatives not implemented for volumes');
        third_der = [];
        fourth_der = [];
      end
    end

  elseif (size(nurbs.knots,2) == 2)
  % NURBS structure represents a surface
    temp = cw(ones(3,1),:,:);
    pnt = cp./temp;
  
    [cup,cuw] = nrbeval (dnurbs{1}, tt);
    tempu = cuw(ones(3,1),:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp,cvw] = nrbeval (dnurbs{2}, tt);
    tempv = cvw(ones(3,1),:,:);
    jac{2} = (cvp-tempv.*pnt)./temp;

% second derivatives
    if (nargout >= 3) 
      if (exist ('dnurbs2'))
        [cuup, cuuw] = nrbeval (dnurbs2{1,1}, tt);
        tempuu = cuuw(ones(3,1),:,:);
        hess{1,1} = (cuup - (2*cup.*tempu + cp.*tempuu)./temp + 2*cp.*tempu.^2./temp.^2)./temp;

        [cvvp, cvvw] = nrbeval (dnurbs2{2,2}, tt);
        tempvv = cvvw(ones(3,1),:,:);
        hess{2,2} = (cvvp - (2*cvp.*tempv + cp.*tempvv)./temp + 2*cp.*tempv.^2./temp.^2)./temp;

        [cuvp, cuvw] = nrbeval (dnurbs2{1,2}, tt);
        tempuv = cuvw(ones(3,1),:,:);
        hess{1,2} = (cuvp - (cup.*tempv + cvp.*tempu + cp.*tempuv)./temp + 2*cp.*tempu.*tempv./temp.^2)./temp;
        hess{2,1} = hess{1,2};
      else
        warning ('nrbdeval: dnurbs2 missing. The second derivative is not computed');
        hess = [];
      end
    end
    if (nargout >= 4) 
      if (exist ('dnurbs3'))
        [cuuup, cuuuw] = nrbeval (dnurbs3{1,1,1}, tt);
        tempuuu = cuuuw(ones(3,1),:,:);
        third_der{1,1,1} = (cuuup - 3*jac{1}.*tempuu - cp.*tempuuu - 3*hess{1,1}.*tempu)./temp;

        [cvvvp, cvvvw] = nrbeval (dnurbs3{2,2,2}, tt);
        tempvvv = cvvvw(ones(3,1),:,:);
        third_der{2,2,2} = (cvvvp - 3*jac{2}.*tempvv - cp.*tempvvv - 3*hess{2,2}.*tempv)./temp;
        
        [cuuvp, cuuvw] = nrbeval (dnurbs3{1,1,2}, tt);
        tempuuv = cuuvw(ones(3,1),:,:);
        third_der{1,1,2} = (cuuvp - 2*jac{1}.*tempuv - jac{2}.*tempuu - cp.*tempuuv - 2*hess{1,2}.*tempu - hess{1,1}.*tempv)./temp;
        third_der{1,2,1} = third_der{1,1,2};
        third_der{2,1,1} = third_der{1,1,2};
  
        [cuvvp, cuvvw] = nrbeval (dnurbs3{1,2,2}, tt);
        tempuvv = cuvvw(ones(3,1),:,:);
        third_der{1,2,2} = (cuvvp - 2*jac{2}.*tempuv - jac{1}.*tempvv - cp.*tempuvv - 2*hess{1,2}.*tempv - hess{2,2}.*tempu)./temp;
        third_der{2,2,1} = third_der{1,2,2};
        third_der{2,1,2} = third_der{1,2,2};
        
      else
        warning ('nrbdeval: dnurbs3 missing. The third derivative is not computed');
        third_der = [];
      end
    end

    if (nargout == 5) 
      if (exist ('dnurbs4'))
        [cuuuup, cuuuuw] = nrbeval (dnurbs4{1,1,1,1}, tt);
        tempuuuu = cuuuuw(ones(3,1),:,:);
        fourth_der{1,1,1,1} = (cuuuup - cp.*tempuuuu - 4*jac{1}.*tempuuu -6*hess{1,1}.*tempuu - 4*third_der{1,1,1}.*tempu)./temp;

        [cvvvvp, cvvvvw] = nrbeval (dnurbs4{2,2,2,2}, tt);
        tempvvvv = cvvvvw(ones(3,1),:,:);
        fourth_der{2,2,2,2} = (cvvvvp - cp.*tempvvvv - 4*jac{2}.*tempvvv -6*hess{2,2}.*tempvv - 4*third_der{2,2,2}.*tempv)./temp;

        [cuuuvp, cuuuvw] = nrbeval (dnurbs4{1,1,1,2}, tt);
        tempuuuv = cuuuvw(ones(3,1),:,:);
        fourth_der{1,1,1,2} = (cuuuvp - cp.*tempuuuv - 3*jac{1}.*tempuuv - jac{2}.*tempuuu -3*hess{1,2}.*tempuu -3*hess{1,1}.*tempuv ...
                                - third_der{1,1,1}.*tempv - 3*third_der{1,1,2}.*tempu)./temp;
                            
        fourth_der{1,1,2,1} = fourth_der{1,1,1,2};
        fourth_der{1,2,1,1} = fourth_der{1,1,1,2};
        fourth_der{2,1,1,1} = fourth_der{1,1,1,2};

        [cuvvvp, cuvvvw] = nrbeval (dnurbs4{1,2,2,2}, tt);
        tempuvvv = cuvvvw(ones(3,1),:,:);
        fourth_der{1,2,2,2} = (cuvvvp - cp.*tempuvvv - 3*jac{2}.*tempuvv - jac{1}.*tempvvv -3*hess{1,2}.*tempvv -3*hess{2,2}.*tempuv ...
                                - third_der{2,2,2}.*tempu - 3*third_der{1,2,2}.*tempv)./temp;
                            
        fourth_der{2,2,1,2} = fourth_der{1,2,2,2};
        fourth_der{2,1,2,2} = fourth_der{1,2,2,2};
        fourth_der{2,2,2,1} = fourth_der{1,2,2,2};

        [cuuvvp, cuuvvw] = nrbeval (dnurbs4{1,1,2,2}, tt);
        tempuuvv = cuuvvw(ones(3,1),:,:);
        fourth_der{1,1,2,2} = (cuuvvp - cp.*tempuuvv - 2*jac{1}.*tempuvv - 2*jac{2}.*tempuuv -hess{1,1}.*tempvv -hess{2,2}.*tempuu ...
                                -4*hess{1,2}.*tempuv - 2*third_der{1,1,2}.*tempv - 2*third_der{1,2,2}.*tempu)./temp;
                            
        fourth_der{1,2,2,1} = fourth_der{1,1,2,2};
        fourth_der{2,2,1,1} = fourth_der{1,1,2,2};
        fourth_der{1,2,1,2} = fourth_der{1,1,2,2};
        fourth_der{2,1,2,1} = fourth_der{1,1,2,2};
        fourth_der{2,1,1,2} = fourth_der{1,1,2,2};
        
        clear cuup cuuw tempuu cvvp cvvw tempvv cuvp cuvw tempuv
        clear cuuup cuuuw tempuuu cvvvp cvvvw tempvvv cuuvp cuuvw tempuuv cuvvp cuvvw tempuvv
        clear cuuuup cuuuuw tempuuuu cvvvvp cvvvvw tempvvvv cuuuvp cuuuvw tempuuuv cuuvvp cuuvvw tempuuvv cvvvup cvvvuw tempvvvu
        
      else
        warning ('nrbdeval: dnurbs4 missing. The fourth derivative is not computed');
        fourth_der = [];
      end
    end

    
  end
else

  % NURBS is a curve
  temp = cw(ones(3,1),:);
  pnt = cp./temp;
  
  % first derivative
  [cup,cuw] = nrbeval (dnurbs,tt);
  temp1 = cuw(ones(3,1),:);
  jac = (cup-temp1.*pnt)./temp;

  % second derivative
  if (nargout >= 3 && exist ('dnurbs2'))
    [cuup,cuuw] = nrbeval (dnurbs2, tt);
    temp2 = cuuw(ones(3,1),:);
    hess = (cuup - (2*cup.*temp1 + cp.*temp2)./temp + 2*cp.*temp1.^2./temp.^2)./temp; 
    if (nargout >= 4 && exist ('dnurbs3')) 
    
      [cuuup, cuuuw] = nrbeval (dnurbs3, tt);
      temp3 = cuuuw(ones(3,1),:);
      third_der = (cuuup - 3*jac.*temp2 - cp.*temp3 - 3*hess.*temp1)./temp;
      if (nargout >= 5 && exist ('dnurbs4')) 
    
        [cuuuup, cuuuuw] = nrbeval (dnurbs4, tt);
        temp4 = cuuuuw(ones(3,1),:);
        fourth_der = (cuuuup - cp.*temp4 - 4*jac.*temp3 -6*hess.*temp2 - 4*third_der.*temp1)./temp;
        if (iscell (tt))
          fourth_der = {fourth_der};
        end

      end
      
      if (iscell (tt))
        third_der = {third_der};
      end      
    end
    
    if (iscell (tt))
      hess = {hess};
    end    
  end 
  
  if (iscell (tt))
    jac = {jac};
  end
  
end

varargout{1} = pnt;
varargout{2} = jac;
if (nargout >= 3)
  varargout{3} = hess;
end
if (nargout >= 4)
  varargout{4} = third_der;
end
if (nargout == 5)
  varargout{5} = fourth_der;
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
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 .5 1 1 1];
%! knots{3} = [0 0 0 0 1 1 1 1];
%! cx = [0 0.5 1]; nx = length(cx);
%! cy = [0 0.25 0.75 1]; ny = length(cy);
%! cz = [0 1/3 2/3 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! ders = nrbderiv(nurbs);
%! [points,jac] = nrbdeval(nurbs,ders,tt);
%! assert(points,tt,1e-10)
%! assert(jac{1}(1,:,:),ones(size(jac{1}(1,:,:))),1e-12)
%! assert(jac{2}(2,:,:),ones(size(jac{2}(2,:,:))),1e-12)
%! assert(jac{3}(3,:,:),ones(size(jac{3}(3,:,:))),1e-12)
%! 
%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 0 1 1 1 1];
%! knots{3} = [0 0 0 1 1 1];
%! cx = [0 0 1]; nx = length(cx);
%! cy = [0 0 0 1]; ny = length(cy);
%! cz = [0 0.5 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! coefs = coefs([2 1 3 4],:,:,:);
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! dnurbs = nrbderiv(nurbs);
%! [points, jac] = nrbdeval(nurbs,dnurbs,tt);
%! assert(points,[y.^3 x.^2 z]',1e-10);
%! assert(jac{2}(1,:,:),3*y'.^2,1e-12)
%! assert(jac{1}(2,:,:),2*x',1e-12)
%! assert(jac{3}(3,:,:),ones(size(z')),1e-12)
