function C = basiskntins (deg, kc, kf)

% Compute the coefficient matrix for non-uniform B-splines subdivision.
%
% This represents the B-spline basis given by a coarse knot vector
%  in terms of the B-spline basis of a finer knot vector.
%
% The function is implemented for the univariate case, based on
%  Algorithm A5.4 from 'The NURBS BOOK' pg164.
%
%
% Calling Sequence:
% 
%    S = basiskntins (deg, kc, kf);
%
%    INPUT:
%   
%      deg - degree of the first knot vector
%      kc  - coarse knot vector
%      kf  - fine knot vector
%   
%    OUTPUT:
%   
%      S - The matrix relating the two spaces, of size (deg-nu, deg-nt) 
%           with nu = numel(u)-deg-1, nt = numel(t)-deg-1
%   
%    Copyright (C) 2015, 2016, 2018 Rafael Vazquez
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

nk = numel(kc);
nc = nk - (deg+1);

u = new_knots(kc, kf);
nu = numel(u);
nf = nc + nu;
if (nu == 0)
  C = speye (nf, nc);
  return
else
  C = sparse (nf, nc);
end

ik = zeros(1,nk+nu);
                                                     
n = nc - 1;
r = nu - 1;

m = nc + deg;
a = findspan(n, deg, u(1), kc);
b = findspan(n, deg, u(end), kc);
b = b+1;

C(1:a-deg+1,1:a-deg+1) = speye(a-deg+1);
C(b+nu:nc+nu,b:nc) = speye(nc-b+1);

ik(1:a+1) = kc(1:a+1);
ik(b+deg+nu+1:m+nu+1) = kc(b+deg+1:m+1);

ii = b + deg - 1;
ss = ii + nu;

for jj=r:-1:0
   ind = (a+1):ii;
   ind = ind(u(jj+1)<=kc(ind+1));
   C(ind+ss-ii-deg,:) = 0;
   C(ind+ss-ii-deg,ind-deg) = speye(numel(ind));
   ik(ind+(ss-ii)+1) = kc(ind+1);
   ii = ii - numel(ind);
   ss = ss - numel(ind);

   C(ss-deg,:) = C(ss-deg+1,:);

   for l=1:deg
       ind = ss - deg + l;
       alfa = ik(ss+l+1) - u(jj+1);
       if abs(alfa) == 0
           C(ind,:) = C(ind+1,:);
       else
           alfa = alfa/(ik(ss+l+1) - kc(ii-deg+l+1));
           C(ind,:) = C(ind,:)*alfa + C(ind+1,:)*(1-alfa);
       end
   end

   ik(ss+1) = u(jj+1);
   ss = ss - 1;
end

end

function u = new_knots (kc, kf)
% Find the new knots, with the correct multiplicity
[valc, multc] = unique (kc, 'last');
multc = diff ([0 multc(:)']);
[valf, multf] = unique (kf, 'last');
multf = diff ([0 multf(:)']);

unew = setdiff (kf, kc);
[~,posf] = ismember (unew, valf);
mult_new = multf(posf);

[urep, indc, indf] = intersect (valc, valf);
mult_rep = multf(indf) - multc(indc);
urep = urep(mult_rep>0);
mult_rep = mult_rep(mult_rep>0);

mult = [mult_new mult_rep];
u = [unew, urep];

ind = zeros (numel(kf)-numel(kc), 1);
ind(cumsum([1 mult(:)'])) = 1;
u = sort (u(cumsum(ind(1:end-1))));

end