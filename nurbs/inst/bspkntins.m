function [ic,ik,C] = bspkntins(d,c,k,u)

% BSPKNTINS:  Insert knots into a B-Spline
%
% Calling Sequence:
% 
%   [ic,ik] = bspkntins(d,c,k,u)
%
%  INPUT:
% 
%    d - spline degree             integer
%    c - control points            double  matrix(mc,nc)      
%    k - knot sequence             double  vector(nk) 
%    u - new knots                 double  vector(nu)               
% 
%  OUTPUT:
% 
%    ic - new control points double  matrix(mc,nc+nu) 
%    ik - new knot sequence  double  vector(nk+nu)
% 
%  Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164.
% 
%    Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton, 2010-2016 Rafael Vazquez
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

[mc,nc] = size(c);
u  = sort(u);
nu = numel(u);
nk = numel(k);
                                                     % 
                                                     % int bspkntins(int d, double *c, int mc, int nc, double *k, int nk,
                                                     %               double *u, int nu, double *ic, double *ik)
                                                     % {
                                                     %   int ierr = 0;
                                                     %   int a, b, r, l, i, j, m, n, s, q, ind;
                                                     %   double alfa;
                                                     %
                                                     %   double **ctrl  = vec2mat(c, mc, nc);
ic = zeros(mc,nc+nu);                                %   double **ictrl = vec2mat(ic, mc, nc+nu);
ik = zeros(1,nk+nu);
                                                     %
n = nc - 1;                                          %   n = nc - 1;
r = nu - 1;                                          %   r = nu - 1;
                                                     %
m = n + d + 1;                                       %   m = n + d + 1;
a = findspan(n, d, u(1), k);                         %   a = findspan(n, d, u[0], k);
b = findspan(n, d, u(r+1), k);                       %   b = findspan(n, d, u[r], k);
b = b+1;                                             %   ++b;
                                                     %
                                                     %   for (q = 0; q < mc; q++)  {
ic(:,1:a-d+1) = c(:,1:a-d+1);                        %     for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q];
ic(:,b+nu:nc+nu) = c(:,b:nc);                        %     for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q];
                                                     %   }
                                                     
ik(1:a+1) = k(1:a+1);                                %   for (j = 0; j <= a; j++)   ik[j] = k[j];
ik(b+d+nu+1:m+nu+1) = k(b+d+1:m+1);                      %   for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];
                                                     %
ii = b + d - 1;                                      %   i = b + d - 1;
ss = ii + nu;                                        %   s = b + d + r;

for jj=r:-1:0                                        %   for (j = r; j >= 0; j--) {
   ind = (a+1):ii;                                   %     while (u[j] <= k[i] && i > a) {
   ind = ind(u(jj+1)<=k(ind+1));                     %       for (q = 0; q < mc; q++)
   ic(:,ind+ss-ii-d) = c(:,ind-d);                   %         ictrl[s-d-1][q] = ctrl[i-d-1][q];
   ik(ind+ss-ii+1) = k(ind+1);                       %       ik[s] = k[i];
   ii = ii - numel(ind);                             %       --i;
   ss = ss - numel(ind);                             %       --s;
                                                     %     }
   ic(:,ss-d) = ic(:,ss-d+1);                        %     ictrl[s-d-1][q] = ictrl[s-d][q];

   for l=1:d                                         %     for (l = 1; l <= d; l++)  {
       ind = ss - d + l;                             %       ind = s - d + l;
       alfa = ik(ss+l+1) - u(jj+1);                  %       alfa = ik[s+l] - u[j];
       if abs(alfa) == 0                             %       if (fabs(alfa) == 0.0)    
           ic(:,ind) = ic(:,ind+1);                  %         for (q = 0; q < mc; q++)
                                                     %           ictrl[ind-1][q] = ictrl[ind][q];
       else                                          %       else  {
           alfa = alfa/(ik(ss+l+1) - k(ii-d+l+1));   %         alfa /= (ik[s+l] - k[i-d+l]);
                      
           tmp = (1-alfa) * ic(:,ind+1);             %         for (q = 0; q < mc; q++)
           ic(:,ind) = alfa*ic(:,ind) + tmp;         %           ictrl[ind-1][q] = alfa*ictrl[ind-1][q]+(1.0-alfa)*ictrl[ind][q];
       end                                           %       }
   end                                               %     }
   %
   ik(ss+1) = u(jj+1);                               %     ik[s] = u[j];
   ss = ss - 1;
end                                                  %   }
                                                     %
                                                     %   freevec2mat(ctrl);
                                                     %   freevec2mat(ictrl);
                                                     %
                                                     %   return ierr;
end                                                  % }

