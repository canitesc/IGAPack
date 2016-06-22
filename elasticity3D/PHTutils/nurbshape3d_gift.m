function [N, dN]=nurbshape3d_gift(u_hat, v_hat, w_hat, wgts, Ce, p, q, r)

%calculate the shape function and first derivatives
%INPUT:
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       v_hat - evaluation point in v-direction reference coordinates (from -1 to 1)
%       w_hat - evaluation point in w-direction reference coordinates (from
%       -1 to 1)
%
%       Ce - local Bezier extraction operator
%       p,q,r - polynomial degrees
%OUTPUT: N - value of the (p+1)*(q+1)*(r+1) shape functions at the evaluation point
%        dN - derivatives of the (p+1)*(q+1)*(r+1) shape functions in physical
%        coordinates (with respect to x)
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives
[B_u, dB_u] = bernstein_basis(u_hat,p);
[B_v, dB_v] = bernstein_basis(v_hat,q);
[B_w, dB_w] = bernstein_basis(w_hat,r);

B = zeros(1, (p+1)*(q+1)*(r+1));
dBxi = zeros(1, (p+1)*(q+1)*(r+1));
dBeta = zeros(1, (p+1)*(q+1)*(r+1));
dBzeta = zeros(1, (p+1)*(q+1)*(r+1));

basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            
            %  Bernstein basis functions;
            B(basisCounter) = B_u(i)*B_v(j)*B_w(k);
            dBxi(basisCounter) = dB_u(i)*B_v(j)*B_w(k);
            dBeta(basisCounter) = B_u(i)*dB_v(j)*B_w(k);
            dBzeta(basisCounter) = B_u(i)*B_v(j)*dB_w(k);
        end
    end
end

%form B-spline basis functions using the Bezier extraction operator
B = B*Ce';
dBxi = dBxi*Ce';
dBeta = dBeta*Ce';
dBzeta = dBzeta*Ce';

N = zeros(1, (p+1)*(q+1)*(r+1));
dN = zeros(3, (p+1)*(q+1)*(r+1));


% Multiply each B-spline function with corresponding weight
N(1,:) = B .* wgts';
dN(1,:) = dBxi .* wgts';
dN(2,:) = dBeta .* wgts';
dN(3,:) = dBzeta .* wgts';

% Compute the sums of B-spline functions
w_sum = sum(N(1,:));
dw_xi = sum(dN(1,:));
dw_eta = sum(dN(2,:));
dw_zeta = sum(dN(3,:));

% Compute NURBS basis functions and its first derivatives in
% local coordinates

dN(1,:) = dN(1,:)/w_sum - N*dw_xi/w_sum^2;
dN(2,:) = dN(2,:)/w_sum - N*dw_eta/w_sum^2;
dN(3,:) = dN(3,:)/w_sum - N*dw_zeta/w_sum^2;

N = N/w_sum;
