function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuLSBracketMP(stiff, rhs, PHTelem, p, q, bound_disp)
%impose Dirichlet boundary conditions for elastic rectangle
%fixed (homogeneous) boundary conditions on the left side
%Neumann (traction) boundary conditions on the right side

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];

%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann_left = [];
neumann_right = [];
neumann_up = [];
neumann_down =[];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);


for i=1:length(PHTelem{16})
    if isempty(PHTelem{16}(i).children)        
        if isempty(PHTelem{16}(i).neighbor_up)
            bcdof_left = [bcdof_left, PHTelem{16}(i).nodesGlobal(up_nodes)];  
        end
     
    end
end

%set the boundary degree of freedom and elements from the 4th patch
for i=1:length(PHTelem{4})
    if isempty(PHTelem{4}(i).children)        
        if isempty(PHTelem{4}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{4}(i).nodesGlobal(right_nodes)];  
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

bcval_right_x = zeros(size(bcdof_right_x));
bcval_right_y = -bound_disp*ones(size(bcdof_right_y));

bcdof = [bcdof_left_x, bcdof_left_y, bcdof_right_y];
bcval = [bcval_left_x, bcval_left_y, bcval_right_y];

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
