function [u, convergence] = nrbinverse (nrb, x, varargin)
%
% NRBINVERSE: compute parametric point starting from physical point by
% inverting the NURBS map with a Newton scheme
%
% Calling Sequence:
%
%   u = nrbinverse (nrb, x)
%   u = nrbinverse (nrb, x, options)
%
%    INPUT:
%
%      nrb     - NURBS object
%      x       - physical point
%      options - options in the FIELD/VALUE format. Possible choices:
%        'u0'      : starting point in the parametric domain for Newton 
%                    (Default = .5 * ones (ndim, 1))
%        'MaxIter' : maximum number of Newton iterations (Default = 10)
%        'Display' : if true the some info are shown (Default = true)
%        'TolX'    : tolerance for the step size in Newton iterations
%                    (Default = 1e-8)
%        'TolFun'  : tolerance for the residual in Newton iterations
%                    (Default = 1e-8)
%
%    OUTPUT:
%
%      u            - the parametric points corresponding to x
%      convergence  - false if the method reached the maximum number
%                     of iteration without converging, true otherwise
%
% Copyright (C) 2016 Jacopo Corno
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

  ndim = numel (nrb.number);
%  % Default options
%  persistent p;
%  p = inputParser ();
%  p.addParameter ('u0', .5*ones(ndim, 1), @(x) validateattributes (x, {'numeric'}, {'numel', ndim, '>=', 0, '<=', 1}));
%  p.addParameter ('MaxIter', 10, @(x) validateattributes (x, {'numeric'}, {'scalar'}));
%  p.addParameter ('Display', true, @(x) validateattributes (x, {'logical'}, {}));
%  p.addParameter ('TolX', 1e-8, @(x) validateattributes (x, {'numeric'}, {'scalar'}));
%  p.addParameter ('TolFun', 1e-8, @(x) validateattributes (x, {'numeric'}, {'scalar'}));
%  p.parse (varargin{:});
%  options = p.Results;
  
  % Default options
  options = struct ('u0'      , .5*ones (ndim, 1), ...
                    'MaxIter' , 10, ...
                    'Display' , true, ...
                    'TolX',     1e-8, ...
                    'TolFun',   1e-8);

  % Read the acceptable names
  optionNames = fieldnames (options);

  % Count arguments
  nargin = length (varargin);
  if (round (nargin/2) ~= nargin/2)
     error ('NRBINVERSE needs propertyName/propertyValue pairs');
  end
  
  % Check options passed
  for pair = reshape (varargin, 2, [])
    if any (strcmp (pair{1}, optionNames))
      options.(pair{1}) = pair{2};
    else
      error('%s is not a recognized parameter name', pair{1});
    end
  end
  
  % x as column vector
  x = x(:);
  
  % Define functions for Newton iteration
  f = @(U) nrbeval (nrb, num2cell (U)) - x;
  jac = @(U) nrbjacobian (nrb, num2cell (U));
  
  % Newton cycle
  u_old = options.u0(:);

  if (iscell (nrb.knots))
    first_knot = reshape (cellfun (@(x) x(1),nrb.knots), size(u_old));
    last_knot = reshape (cellfun (@(x) x(end),nrb.knots), size(u_old));
  else
    first_knot = nrb.knots(1);
    last_knot = nrb.knots(end);
  end
  convergence = false;
  
  for iter = 1:options.MaxIter

    u_new = u_old - jac (u_old) \ f (u_old);

    % Check if the point is outside the parametric domain
    u_new = max (u_new, first_knot);
    u_new = min (u_new, last_knot);
    
    % Error control
    if (norm (u_new - u_old) < options.TolX && norm (f (u_new)) < options.TolFun)
      if (options.Display)
        fprintf ('Newton scheme converged in %i iteration.\n', iter);
      end
      convergence = true;
      break;
    end
    
    u_old = u_new;
    
  end

  if (~convergence)
    fprintf ('Newton scheme reached the maximum number of iterations (%i) without converging.\n', options.MaxIter);
  end
  u = u_new;

end

function jac = nrbjacobian (nrb, u)

  ders = nrbderiv (nrb);
  [~, jac] = nrbdeval (nrb, ders, u);
  jac = [jac{:}];  

end

%!test
%! nrb = nrb4surf ([0 0], [1 0], [2 3], [5 4]);
%! p = nrbeval (nrb, {.25 .75});
%! u = nrbinverse (nrb, p, 'Display', false);
%! assert (norm (u - [.25; .75]) < 1e-8);
%!
%!test
%! nrb = nrb4surf ([0 0], [1 0], [2 3], [5 4]);
%! nrb = nrbdegelev (nrbextrude (nrb, [0 2 1]), [3 3 3]);
%! p = nrbeval (nrb, {.25 .75 .05});
%! u = nrbinverse (nrb, p, 'Display', false, 'TolX', 1e-12, 'TolFun', 1e-10);
%! assert (norm (u - [.25; .75; .05]) < 1e-8);
%!
