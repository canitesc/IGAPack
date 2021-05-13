function mnrb = nrbmodp (nrb, move, index)

%
% NRBMODP: Modify the coordinates of specific control points of any NURBS
% map. The weight is not changed.
%
% Calling Sequence:
% 
%   nrb = nrbmodp (nrb, move, index);
%   
%    INPUT:
%   
%      nrb   - NURBS map to be modified.
%      move  - vector specifying the displacement of all the ctrl points.
%      index - indeces of the control points to be modified.
%   
%    OUTPUT:
%   
%      mnrb - the modified NURBS.
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
  
  move = reshape (move, 3, 1);

  mnrb = nrb;
  [ii, jj, kk] = ind2sub (nrb.number, index);
  for count = 1:numel (ii)
    mnrb.coefs(1:3,ii(count),jj(count),kk(count)) = nrb.coefs(1:3,ii(count),jj(count),kk(count)) + ...
      move * nrb.coefs(4,ii(count),jj(count),kk(count));
  end
  
end

