function [ stiff, rhs, bcdof, bcval ] = imposeDirichletCantileverBar(stiff, rhs, PHTelem, p, q, r, bound_disp)
%impose Dirichlet boundary conditions for the cantilever bar
%fixed (homogeneous) boundary conditions on the bottom
%displacement (bound_disp units) in the up direction on the end of the vertical bar

%detect which nodes have support on the left and right boundary
bcdof_down = [];
bcdof_up = [];

%numPatches = length(PHTelem);

%define side node indices
down_nodes = 1:(p+1)*(q+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);

for bottomPatchIndex=1:9
    for i=1:length(PHTelem{bottomPatchIndex})
        if isempty(PHTelem{bottomPatchIndex}(i).children)
            if isempty(PHTelem{bottomPatchIndex}(i).neighbor_down)
                bcdof_down = [bcdof_down, PHTelem{bottomPatchIndex}(i).nodesGlobal(down_nodes)];
            end
        end
    end
end

for topPatchIndex=10%[1:4,6:10]
    for i=1:length(PHTelem{topPatchIndex})
        if isempty(PHTelem{topPatchIndex}(i).children)
            if isempty(PHTelem{topPatchIndex}(i).neighbor_up)
                bcdof_up = [bcdof_up, PHTelem{topPatchIndex}(i).nodesGlobal(up_nodes)];
            end
        end
    end
end

%remove duplicated entries
bcdof_down = unique(bcdof_down);
bcdof_up = unique(bcdof_up);

%take into account that there are 3 global indices (displacements) per node
bcdof_down_x = 3*bcdof_down-2;
bcdof_down_y = 3*bcdof_down-1;
bcdof_down_z = 3*bcdof_down;

bcdof_up_x = 3*bcdof_up-2;
bcdof_up_y = 3*bcdof_up-1;
bcdof_up_z = 3*bcdof_up;

%impose the prescribed displacement for each global index
bcval_down_x = zeros(size(bcdof_down_x));
bcval_down_y = zeros(size(bcdof_down_y));
bcval_down_z = zeros(size(bcdof_down_z));


bcval_up_x = zeros(size(bcdof_up_x));
bcval_up_y = bound_disp*ones(size(bcdof_up_y));
bcval_up_z = zeros(size(bcdof_up_z));


% bcdof = [bcdof_down_x, bcdof_down_y, bcdof_down_z, bcdof_up_x, bcdof_up_y, bcdof_up_z];
% bcval = [bcval_down_x, bcval_down_y, bcval_down_z, bcval_up_x, bcval_up_y, bcval_up_z];

bcdof = [bcdof_down_x, bcdof_down_y, bcdof_down_z, bcdof_up_y];
bcval = [bcval_down_x, bcval_down_y, bcval_down_z, bcval_up_y];

[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);