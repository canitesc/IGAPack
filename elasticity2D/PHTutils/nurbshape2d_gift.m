function [N, dN]=nurbshape2d_gift(u_hat, v_hat, wgts, Ce, p, q)

%calculate the shape function and first derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_h - evaluation point in v-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p,q - polynomial degrees
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the p+1 shape functions in physical
%        coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives 
[B_u, dB_u] = bernstein_basis(u_hat,p);
[B_v, dB_v] = bernstein_basis(v_hat,q);

B = zeros(1, (p+1)*(q+1));
dBxi = zeros(1, (p+1)*(q+1));
dBeta = zeros(1, (p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        
        %  Bernstein basis functions;    
        B(basisCounter) = B_u(i)*B_v(j);
        dBxi(basisCounter) = dB_u(i)*B_v(j);
        dBeta(basisCounter) = B_u(i)*dB_v(j);        
    end
end

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';
dBeta = dBeta*Ce';

N = zeros(1, (p+1)*(q+1));
dN = zeros(2, (p+1)*(q+1));


% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));

% Compute NURBS basis functions and its first derivatives in
% local coordinates

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;

N = N/w_sum;
