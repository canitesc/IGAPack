function [N, coord, dNlocal, dN, detj]=basisFun(u_hat,cpts,w,Ce,p)

%calculate the shape function and second derivatives
%       u_hat - evaluation point in reference coordinates (from -1 to 1)
%       cpts - control points
%       w - weight
%       Ce - local Bezier extraction operator
%       p - polynomial degree
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        coord - the coordinates of the evaluation point in physical space
%        dNlocal - derivatives of the p+1 shape functions in local coordinates
%        (with respect to xi)
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
%       detj - determinant of the jacobian of the transfromation (dx/dxi)
%       ddN - 2nd derivatives of the p+1 shape functions in physical
%       coordinates (with respect to x)

%nsd = 1; 
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[B, dB, ddB] = bernstein_basis(u_hat,p);
M=[B;dB;ddB]*Ce'; 

B = M(1,:);
dBxi = M(2,:);
%ddBxi = M(3,:);

% Multiply each B-spline function with corresponding weight
N = B .* w';
dN = dBxi .* w';
%ddN = ddBxi .* w';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
%d2w_xi = sum(ddN(1,:));


% Compute NURBS basis functions and its first and second derivatives in
% local coordinates
%ddN = ddN/w_sum - (2*dw_xi*dN + N*d2w_xi)/w_sum^2 + 2*N*dw_xi^2/w_sum^3;
dN = dN/w_sum - N*dw_xi/w_sum^2;
N = N/w_sum;

coord = N*cpts;

% Compute Jacobian matrix
dxdxi = dN*cpts;

% Set up the Hessian and the matrix of squared first derivatives
%d2xdxi2 = ddN*cpts;
%dxdxi2 = dxdxi^2;

% Solve for first derivatives in global coordinates 
dNlocal = dN;
dN = dxdxi\dN;

% Solve for second derivatives in global coordinates 
%ddN = dxdxi2\(ddN - d2xdxi2*dN);
detj = det(dxdxi);

