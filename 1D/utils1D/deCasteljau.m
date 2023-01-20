function [ Ce1, Ce2, newB ] = deCasteljau( Ce, knot1, knot2, cpts, w)
% Split element with Bezier extraction operator Ce into 2 elements with
% Bezier extraction operator Ce1, Ce2, calculates new control points and
% weights, solves a linear system with 2p-4 unknows on each patch

num_rows = size(Ce,1);
p = num_rows - 1;


Ce1 = zeros(num_rows, num_rows);
Ce2 = zeros(num_rows, num_rows);

temp = zeros(num_rows, num_rows);
for i=1:size(Ce,1)
    cur_row = Ce(i,:);
    temp(1,:) = cur_row;
    for j=1:num_rows-1
        temp(j+1,1:end-1) = (temp(j,1:end-1)+temp(j,2:end))/2;
    end
    Ce1(i,:) = temp(:,1);
    for j=1:num_rows
        Ce2(i,j) = temp(end+1-j,j);
    end
    
    %zero out entries corresponding to the new knot   
    Ce1(:, end-1:end) = 0;
    Ce2(:, 1:2) = 0;      
end

knot_mid = (knot1+knot2)/2;

newKnotVector = [zeros(1,p+1), knot1*ones(1,p-1), knot_mid*ones(1,p-1), knot2*ones(1,p-1), ones(1, p+1)];
[newC, ~] = bezierExtraction(newKnotVector,p);



%comment this out to get different basis functions
Ce1(end-1:end, :) = newC(end-1:end,:,2);
Ce2(1:2, :) = newC(1:2,:,3);

%vs.

%Ce1(end-p+2:end, :) = newC(end-p+2:end,:,2);
%Ce2(1:p-1, :) = newC(1:p-1,:,3);

%calculate the new control points        
%pick the sample points for interpolating

[~,xi_patch] = quadrature(2*p-4,'GAUSS',1);
xi_patch = sort(xi_patch);
int_lhs = zeros(2*p-4,2*p-2);
int_rhs = zeros(2*p-4,2);

%calculate the matrix system rhs
for i=1:2*p-4
    u_hat = xi_patch(i);     
    [N, ~, N_proj, ~]=nurbshape1d_proj(u_hat,cpts,w,Ce,p);
    int_rhs(i,1) = N_proj;
    int_rhs(i,2) = N;
    
    if u_hat <= 0
        %map the point from [-1, 0] to [-1,1]
        u_hat_ref = 2*u_hat + 1;
        B = bernstein_basis(u_hat_ref,p);
        M = B*Ce1';
        int_lhs(i,1:p+1) = M;
    else
        %map the point from (0, 1] to (-1, 1]
        u_hat_ref = 2*u_hat-1;
        B = bernstein_basis(u_hat_ref,p);
        M = B*Ce2';
        int_lhs(i,p:2*p) = M;
    end    
end

red_vector = [cpts(1:2).*w(1:2), w(1:2); zeros(2*p-4,2); cpts(end-1:end).*w(end-1:end), w(end-1:end)];
int_rhs = int_rhs - int_lhs*red_vector;

int_lhs(:,[1:2,end-1:end]) = [];
newB = int_lhs\int_rhs;
%divide through by the weights
newB(:,1) = newB(:,1)./newB(:,2);