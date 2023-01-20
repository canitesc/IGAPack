function [N, dNlocal, N_proj, dNlocal_proj]=nurbshape1d_proj(u_hat,cpts,w,Ce,p)

%calculate the (projective) shape function in parameter space
%INPUT: 
%       u_hat - evaluation point in reference coordinates (from -1 to 1)
%      
%       Ce - local Bezier extraction operator
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        coord - the coordinates of the evaluation point in physical space
%        dNlocal - derivatives of the p+1 shape functions in local coordinates
%        (with respect to xi)
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
%       detj - determinant of the jacobian of the transfromation (dx/dxi)
%       ddN - 2nd derivatives of the p+1 shape functions in physical
%       coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[B, dB] = bernstein_basis(u_hat,p);
M=[B;dB]*Ce'; 

B = M(1,:);
dBxi = M(2,:);

% Multiply each B-spline function with corresponding weight
N = B .* w';
dNlocal = dBxi .* w';

N_proj = N*cpts;
dNlocal_proj = dNlocal*cpts;

N = sum(N);
dNlocal = sum(dNlocal);