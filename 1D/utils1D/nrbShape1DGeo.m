function [N, dN]=nrbShape1DGeo(u_hat, wgts, Ce, p)

%calculate the shape function and first derivatives
%INPUT:
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)%
%       Ce - local Bezier extraction operator
%       p - polynomial degree
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives
[B, dBxi] = bernstein_basis(u_hat,p);

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';

% Multiply each B-spline function with corresponding weight
N = B .* wgts';
dN = dBxi .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N);
dw_xi = sum(dN);

% Compute NURBS basis functions and its first derivatives in
% local coordinates
dN = dN/w_sum - N*dw_xi/w_sum^2;
N = N/w_sum;
