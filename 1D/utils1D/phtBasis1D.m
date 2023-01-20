function [Bu, dBu] = phtBasis1D( u_hat, Ce, p )

%calculate the shape function and first derivatives
%INPUT: 
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p - polynomial degree
%OUTPUT: Bu - value of the p+1 shape functions at the evaluation point
%        dBu - derivatives of the shape functions in %       
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[Bu, dBu] = bernstein_basis(u_hat,p);

%form B-spline basis functions using the Bezier extraction operator
Bu = Bu*Ce';
dBu = dBu*Ce';
