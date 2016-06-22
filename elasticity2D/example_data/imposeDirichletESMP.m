function [ stiff, rhs, bcdof, bcval ] = imposeDirichletESMP(stiff, rhs, PHTelem, p, q, bound_disp)
%impose Dirichlet boundary conditions for elastic square
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side
%supports multipatches. Assumes boundary conditions are imposed on the
%first and last patches

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];

%define side node indices
%down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
%up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

numPatches = length(PHTelem);

%determine the bcdofs on the fixed edge (left edge of the first patch)
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{1}(i).nodesGlobal(left_nodes)];
        end
    end
end

%determine the bcdofs on the displaced edge (right edge of the last patch)
for i=1:length(PHTelem{numPatches})
    if isempty(PHTelem{numPatches}(i).children)
        if isempty(PHTelem{numPatches}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{numPatches}(i).nodesGlobal(right_nodes)];
        end
    end
end


%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

%take into account that there are 2 global indices (displacements) per node
bcdof_right_x = 2*bcdof_right-1;
bcdof_right_y = 2*bcdof_right;

bcdof_left_x = 2*bcdof_left-1;
bcdof_left_y = 2*bcdof_left;

%impose the prescribed displacement for each global index
bcval_left_x = zeros(size(bcdof_left_x));
bcval_left_y = zeros(size(bcdof_left_y));

bcval_right_x = bound_disp*ones(size(bcdof_right_x));
bcval_right_y = zeros(size(bcdof_right_y));

bcdof = [bcdof_left_x, bcdof_left_y, bcdof_right_x, bcdof_right_y];
bcval = [bcval_left_x, bcval_left_y, bcval_right_x, bcval_right_y];

[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
