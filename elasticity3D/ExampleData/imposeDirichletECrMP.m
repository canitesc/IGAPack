function [ stiff, rhs, bcdof, bcval ] = imposeDirichletECrMP(stiff, rhs, PHTelem, p, q, r, bound_disp)
%impose Dirichlet boundary conditions for elastic crack
%fixed (homogeneous) boundary conditions on the bottom side
%displacement (bound_disp units) to the top on the top side

%detect which nodes have support on the left and right boundary
bcdof_bottom = [];
bcdof_top = [];

numPatches = length(PHTelem);

%define side node indices
down_nodes = 1:(p+1)*(q+1);
%right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
%left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_down) && (PHTelem{1}(i).vertex(3)==0)
            bcdof_bottom = [bcdof_bottom, PHTelem{1}(i).nodesGlobal(down_nodes)];
        end
    end
end

for i=1:length(PHTelem{2})
    if isempty(PHTelem{2}(i).children)
        if isempty(PHTelem{2}(i).neighbor_down) && (PHTelem{2}(i).vertex(3)==0)
            bcdof_bottom = [bcdof_bottom, PHTelem{2}(i).nodesGlobal(down_nodes)];
        end
    end
end

for i=1:length(PHTelem{3})
    if isempty(PHTelem{3}(i).children)
        if isempty(PHTelem{3}(i).neighbor_up) && (PHTelem{3}(i).vertex(6)==1)
            bcdof_top = [bcdof_top, PHTelem{3}(i).nodesGlobal(up_nodes)];
        end
    end
end

for i=1:length(PHTelem{4})
    if isempty(PHTelem{4}(i).children)
        if isempty(PHTelem{4}(i).neighbor_up) && (PHTelem{4}(i).vertex(6)==1)
            bcdof_top = [bcdof_top, PHTelem{4}(i).nodesGlobal(up_nodes)];
        end
    end
end



%remove duplicated entries
bcdof_bottom = unique(bcdof_bottom);
bcdof_top = unique(bcdof_top);

%take into account that there are 3 global indices (displacements) per node
bcdof_top_x = 3*bcdof_top-2;
bcdof_top_y = 3*bcdof_top-1;
bcdof_top_z = 3*bcdof_top;

bcdof_bottom_x = 3*bcdof_bottom-2;
bcdof_bottom_y = 3*bcdof_bottom-1;
bcdof_bottom_z = 3*bcdof_bottom;

%impose the prescribed displacement for each global index
bcval_bottom_x = zeros(size(bcdof_bottom_x));
bcval_bottom_y = zeros(size(bcdof_bottom_y));
bcval_bottom_z = zeros(size(bcdof_bottom_z));


bcval_top_x = zeros(size(bcdof_top_x));
bcval_top_y = zeros(size(bcdof_top_y));
bcval_top_z = bound_disp*ones(size(bcdof_top_z));


bcdof = [bcdof_bottom_x, bcdof_bottom_y, bcdof_bottom_z, bcdof_top_x, bcdof_top_y, bcdof_top_z];
bcval = [bcval_bottom_x, bcval_bottom_y, bcval_bottom_z, bcval_top_x, bcval_top_y, bcval_top_z];

[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);