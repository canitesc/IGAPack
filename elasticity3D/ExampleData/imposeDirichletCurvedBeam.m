function [ stiff, rhs, bcdof, bcval ] = imposeDirichletCurvedBeam(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_press)
%impose Dirichlet and Neumann boundary conditions for the hollow sphere
%symmetry condition on the flat edges, pressure on the interior surface


%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_front = [];
bcdof_right = [];
bcdof_down = [];
bcdof_back = [];
bcdof_left_back = [];
bcdof_up = [];

%define side node indices
down_nodes = 1:(p+1)*(q+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

front_nodes = [];
neumann_down = [];
neumann_left = [];

for i=1:r+1
    front_nodes = [front_nodes, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end
%front_nodes
back_nodes = front_nodes + (p+1)*q;
back_left_nodes = intersect(back_nodes, left_nodes);

for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        if isempty(PHTelem(i).neighbor_front)
            bcdof_front = [bcdof_front, PHTelem(i).nodes(front_nodes)];
        end
        if isempty(PHTelem(i).neighbor_down)
            neumann_down = [neumann_down; i, 5];
            bcdof_down = [bcdof_down, PHTelem(i).nodes(down_nodes)];
        end
        if isempty(PHTelem(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem(i).nodes(up_nodes)];
        end
        
        if isempty(PHTelem(i).neighbor_left) && isempty(PHTelem(i).neighbor_back)
            bcdof_left_back = [bcdof_left_back, PHTelem(i).nodes(back_left_nodes)];
         
        end
        if isempty(PHTelem(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem(i).nodes(right_nodes)];
        end
        if isempty(PHTelem(i).neighbor_back)
            bcdof_back = [bcdof_back, PHTelem(i).nodes(back_nodes)];
        end
    end
end

%remove duplicated entries
bcdof_front = unique(bcdof_front);
%bcdof_front = setdiff(bcdof_front,bcdof_left);

bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);
bcdof_back = unique(bcdof_back);
%bcdof_back = setdiff(bcdof_back,bcdof_left);

%bcdof_down = unique(bcdof_down);

%take into account that there are 3 global indices (displacements) per node
bcdof_up_x = 3*bcdof_up-2;
bcdof_up_y = 3*bcdof_up-1;
bcdof_up_z = 3*bcdof_up;
%
bcdof_down_x = 3*bcdof_down-2;
bcdof_down_y = 3*bcdof_down-1;
bcdof_down_z = 3*bcdof_down;

bcdof_front_x = 3*bcdof_front-2;
bcdof_front_y = 3*bcdof_front-1;
bcdof_front_z = 3*bcdof_front;

bcdof_back_x = 3*bcdof_back-2;
bcdof_back_y = 3*bcdof_back-1;
bcdof_back_z = 3*bcdof_back;

bcdof_left_back_y = 3*bcdof_left_back-1;

% %impose the symmetry boundary conditions
% bcval_left_x = zeros(size(bcdof_left_x));
% bcval_left_y = zeros(size(bcdof_left_y));
% bcval_left_z = zeros(size(bcdof_left_z));
%
%
% bcval_front_x = zeros(size(bcdof_front_x));
% bcval_front_y = zeros(size(bcdof_front_y));
% bcval_front_z = zeros(size(bcdof_front_z));
%
% bcval_right_x = zeros(size(bcdof_right_x));
% bcval_right_y = zeros(size(bcdof_right_y));
% bcval_right_z = zeros(size(bcdof_right_z));

% size(bcdof_left_y)
% size(bcdof_front_z)
% size(bcdof_right_x)
% size(bcdof_back_x)
% size(bcdof_back_y)
%

bound_disp = 0.1;
bcdof = [bcdof_front_x, bcdof_front_z, bcdof_back_x, bcdof_back_z, bcdof_left_back_y];
bcval = [bound_disp*ones(size(bcdof_front_x)), zeros(size(bcdof_front_z)), zeros(size(bcdof_back_x)), zeros(size(bcdof_back_z)), zeros(size(bcdof_left_back_y))];
%bcval = [bcval_left_y, bcval_front_z, bcval_right_x];
%bcval = zeros(size(bcdof));



[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);

