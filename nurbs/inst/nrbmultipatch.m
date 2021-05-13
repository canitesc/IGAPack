function [interfaces, boundary] = nrbmultipatch (nurbs, tol)

%
% NRBMULTIPATCH: construct the information for gluing conforming NURBS patches, using the same format as in GeoPDEs.
% 
% Calling Sequence:
% 
%   [interfaces, boundary] = nrbmultipatch (nurbs, [tol]);
% 
% INPUT:
% 
%   nurbs   : an array of NURBS surfaces or volumes (not both), see nrbmak.
%   tol     : relative tolerance to compare knots and control points at the interfaces.
% 
% OUTPUT: 
% 
%   interfaces: array with the information for each interface, that is:
%      - number of the first patch (patch1), and the local side number (side1)
%      - number of the second patch (patch2), and the local side number (side2)
%      - flag (faces and volumes), ornt1, ornt2 (only volumes): information
%        on how the two patches match, see below.
%   boundary:   array with the boundary faces that do not belong to any interface
%      - nsides:  total number of sides on the boundary array (numel(boundary))
%      - patches: number of the patch to which the boundary belongs
%      - sides:   number of the local side on the patch
%
% The faces of two patches must match conformingly: the control points must be the same,
%  with the knot vectors (in each direction) related by an affine transformation.
%
% The boundary faces are stored separately, that is, nsides=1 for each boundary.
%  To join several faces under the same condition, the user should do it by hand.
% 
%    Copyright (C) 2014, 2015, 2016, 2017 Rafael Vazquez
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

if (nargin < 2)
  tol = 1e-13;
end

npatch = numel (nurbs);
if (~iscell (nurbs(1).knots))
  ndim = 1;
  compare_sides = @(nrb1, nrb2) max(abs(nrb1.coefs - nrb2.coefs)) < tol;
elseif (size(nurbs(1).knots,2) == 2)
  ndim = 2;
  compare_sides = @(nrb1, nrb2) compare_sides_univariate (nrb1, nrb2, tol);
elseif (size(nurbs(1).knots,2) == 3)
  ndim = 3;
  compare_sides = @(nrb1, nrb2) compare_sides_bivariate (nrb1, nrb2, tol);
end

non_set_faces = cell (npatch, 1);
for ii = 1:npatch
  if (~iscell (nurbs(ii).knots))
    if (ndim ~= 1)
      error ('All the patches must have the same dimension (at least for now)')
    end
  elseif (ndim ~= size(nurbs(ii).knots,2))
    error ('All the patches must have the same dimension (at least for now)')      
  end
  non_set_faces{ii} = 1:2*ndim;
end

num_interfaces = 0;

num_boundaries = 0;
boundary = struct ('nsides', 0, 'patches', [], 'faces', []);

for i1 = 1:npatch
  nrb_faces1 = nrbextract (nurbs(i1));
  for j1 = non_set_faces{i1}
    
% This is to fix a bug when two faces of the same patch form an interface
%  (for instance, in a ring or a torus)
    if (isempty (intersect (non_set_faces{i1}, j1))); continue; end

    nrb1 = nrb_faces1(j1);
%     corners1 = face_corners (nrb1);

    non_set_faces{i1} = setdiff (non_set_faces{i1}, j1);
    flag = 0;

    i2 = i1 - 1;
    while (~flag && i2 < npatch)
      i2 = i2 + 1;
      nrb_faces2 = nrbextract (nurbs(i2));
      j2 = 0;
      while (~flag && j2 < numel (non_set_faces{i2}))
        j2 = j2 + 1;
        nrb2 = nrb_faces2(non_set_faces{i2}(j2));

        if (ndim == 1)
          flag = compare_sides (nrb1, nrb2);
          MsgFlag = false;
        elseif (ndim == 2)
          [flag, MsgFlag] = compare_sides (nrb1, nrb2);
        elseif (ndim == 3)
          [flag, ornt1, ornt2, MsgFlag] = compare_sides (nrb1, nrb2);
        end

        if (MsgFlag)
          display_warning (MsgFlag, i1, j1, i2, j2);
        end
        
      end
    end

    if (flag)
      intrfc.patch1 = i1;
      intrfc.side1 = j1;
      intrfc.patch2 = i2;
      intrfc.side2 = non_set_faces{i2}(j2);
      if (ndim ==3)
        intrfc.flag = flag;
        intrfc.ornt1 = ornt1;
        intrfc.ornt2 = ornt2;
      elseif (ndim == 2)
        intrfc.ornt = flag;
      end

      non_set_faces{i2} = setdiff (non_set_faces{i2}, non_set_faces{i2}(j2));
      num_interfaces = num_interfaces + 1;
      interfaces(num_interfaces) = intrfc;
    else
      bndry.nsides = 1;
      bndry.patches = i1;
      bndry.faces = j1;
      num_boundaries = num_boundaries + 1;
      boundary(num_boundaries) = bndry;
    end
  end
end

if (num_interfaces == 0)
   interfaces = []; 
end
if (num_boundaries == 0)
  boundary = [];
end

end



% Compare the sides of two volumes
function [flag, ornt1, ornt2, MsgFlag] = compare_sides_bivariate (nrb1, nrb2, tol)
  MsgFlag = 0;

  face_corners = @(x) reshape (x.coefs(:, [1 end], [1 end]), 4, []);

  coefs1 = face_corners (nrb1);
  coefs2 = face_corners (nrb2);

% Sort of relative error
  tolcp = tol * max(abs(coefs1(1:3,1) - coefs1(1:3,end)));
  tolknt = tol;

% Should use some sort of relative error
  if (max (max (abs (coefs1 - coefs2))) < tolcp)
    flag = 1; ornt1 = 1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[1 3 2 4])))) < tolcp)
    flag = -1; ornt1 = 1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[3 1 4 2])))) < tolcp)
    flag = -1; ornt1 = -1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[2 1 4 3])))) < tolcp)
    flag = 1; ornt1 = -1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[4 3 2 1])))) < tolcp)
    flag = 1; ornt1 = -1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[4 2 3 1])))) < tolcp)
    flag = -1; ornt1 = -1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[2 4 1 3])))) < tolcp)
    flag = -1; ornt1 = 1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[3 4 1 2])))) < tolcp)
    flag = 1; ornt1 = 1; ornt2 = -1;
  else
    flag = 0; ornt1 = 0; ornt2 = 0;
  end
  
% Reorder control points and knot vectors, to make comparisons easier
  if (flag)
    if (flag == -1)
      nrb2 = nrbtransp (nrb2);
    end
    if (ornt1 == -1)
      nrb2 = nrbreverse (nrb2, 1);
    end
    if (ornt2 == -1)
      nrb2 = nrbreverse (nrb2, 2);
    end
      
    if (nrb1.order ~= nrb2.order)
      flag = 0;
      MsgFlag = -3;
    elseif (nrb1.number ~= nrb2.number)
      flag = 0;
      MsgFlag = -1;
    elseif (any (cellfun (@numel, nrb1.knots) ~= cellfun (@numel, nrb2.knots))) % This is redundant
      flag = 0;
      MsgFlag = -4;
    else
% Pass the knots to the [0 1] interval to compare
      pass_to_01 = @(x) (x - x(1)) / (x(end) - x(1));
      knt1 = cellfun (pass_to_01, nrb1.knots, 'UniformOutput', false);
      knt2 = cellfun (pass_to_01, nrb2.knots, 'UniformOutput', false);
      if (max (abs (nrb1.coefs(:) - nrb2.coefs(:))) > tolcp)
        flag = 0;
        MsgFlag = -2;
      elseif ((max (abs (knt1{1} - knt2{1})) > tolknt) || (max (abs (knt1{2} - knt2{2})) > tolknt))
        flag = 0;
        MsgFlag = -5;
      end      
    end
  end

end



% Compare the sides of two surfaces
function [flag, MsgFlag] = compare_sides_univariate (nrb1, nrb2, tol)
  MsgFlag = 0;
  face_corners = @(x) reshape (x.coefs(:, [1 end]), 4, []);

  coefs1 = face_corners (nrb1);
  coefs2 = face_corners (nrb2);
  
% Sort of relative error
  tolcp = tol * max(abs(coefs1(1:3,1) - coefs1(1:3,end)));
  tolknt = tol;

  if (max (max (abs (coefs1 - coefs2))) < tolcp)
    flag = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[end 1])))) < tolcp)
    flag = -1;
  else
    flag = 0;
  end

  if (flag)
% Reorder control points and knot vectors, to make comparisons easier
    if (flag == -1)
      nrb2 = nrbreverse (nrb2);
    end

    if (nrb1.order ~= nrb2.order)
      flag = 0;
      MsgFlag = -3;
    elseif (nrb1.number ~= nrb2.number)
      flag = 0;
      MsgFlag = -1;
    elseif (numel(nrb1.knots) ~= numel(nrb2.knots)) % This is redundant
      flag = 0;
      MsgFlag = -4;
    else
% Pass the knots to the [0 1] interval to compare
      knt1 = (nrb1.knots - nrb1.knots(1)) / (nrb1.knots(end) - nrb1.knots(1));
      knt2 = (nrb2.knots - nrb2.knots(1)) / (nrb2.knots(end) - nrb2.knots(1));
      if (max (abs (nrb1.coefs(:) - nrb2.coefs(:))) > tolcp)
        flag = 0;
        MsgFlag = -2;
      elseif (max (abs (knt1 - knt2)) > tolknt)
        flag = 0;
        MsgFlag = -5;
      end
    end
  end
  
end

function display_warning (MsgFlag, patch1, face1, patch2, face2)

  switch MsgFlag
    case {-1}
      warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the number ' ... 
                'of control points is different. No information is saved in this case'], patch1, face1, patch2, face2)
    case {-2}
      warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the internal ' ... 
                'control points do not. No information is saved in this case'], patch1, face1, patch2, face2)
    case {-3}
      warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the degree ' ... 
                'is different. No information is saved in this case'], patch1, face1, patch2, face2)
    case {-4}
      warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the number ' ... 
                'of knots is different. No information is saved in this case'], patch1, face1, patch2, face2)
    case {-5}
      warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the ' ... 
                'knot vectors are different. No information is saved in this case'], patch1, face1, patch2, face2)
  end
end