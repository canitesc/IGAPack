function [ stiff, rhs, bcdof, bcval ] = imposeDirichletECUV(stiff, rhs, PHTelem, p, q, r, bound_disp)
%impose Dirichlet boundary conditions for elastic cube
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];

%numPatches = length(PHTelem);

%define side node indices
%down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
%up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

for leftPatchIndex=[1,3]
    for i=1:length(PHTelem{leftPatchIndex})
        if isempty(PHTelem{leftPatchIndex}(i).children)
            if isempty(PHTelem{leftPatchIndex}(i).neighbor_left)
                bcdof_left = [bcdof_left, PHTelem{leftPatchIndex}(i).nodesGlobal(left_nodes)];
            end
        end
    end
end

for rightPatchIndex=[2,4]
    for i=1:length(PHTelem{rightPatchIndex})
        if isempty(PHTelem{rightPatchIndex}(i).children)
            if isempty(PHTelem{rightPatchIndex}(i).neighbor_right)
                bcdof_right = [bcdof_right, PHTelem{rightPatchIndex}(i).nodesGlobal(right_nodes)];
            end
        end
    end
end

%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

%take into account that there are 3 global indices (displacements) per node
bcdof_right_x = 3*bcdof_right-2;
bcdof_right_y = 3*bcdof_right-1;
bcdof_right_z = 3*bcdof_right;

bcdof_left_x = 3*bcdof_left-2;
bcdof_left_y = 3*bcdof_left-1;
bcdof_left_z = 3*bcdof_left;

%impose the prescribed displacement for each global index
bcval_left_x = zeros(size(bcdof_left_x));
bcval_left_y = zeros(size(bcdof_left_y));
bcval_left_z = zeros(size(bcdof_left_z));


bcval_right_x = bound_disp*ones(size(bcdof_right_x));
bcval_right_y = zeros(size(bcdof_right_y));
bcval_right_z = zeros(size(bcdof_right_z));


bcdof = [bcdof_left_x, bcdof_left_y, bcdof_left_z, bcdof_right_x, bcdof_right_y, bcdof_right_z];
bcval = [bcval_left_x, bcval_left_y, bcval_left_z, bcval_right_x, bcval_right_y, bcval_right_z];

[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);