function [ stiff, mass, bcdof ] = imposeDirichletCircPlate(stiff, mass, PHTelem, p, q, r)
%impose Dirichlet boundary conditions for elastic cube
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_front = [];
bcdof_back = [];

numPatches = length(PHTelem);

%define side node indices
%down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
%up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);
front_nodes = [];

for i=1:r+1
    front_nodes = [front_nodes, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end
%front_nodes
back_nodes = front_nodes + (p+1)*q;


for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{1}(i).nodesGlobal(left_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem{numPatches}(i).nodesGlobal(right_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_front)
            bcdof_front = [bcdof_front, PHTelem{numPatches}(i).nodesGlobal(front_nodes)];
        end
        if isempty(PHTelem{1}(i).neighbor_back)
            bcdof_back = [bcdof_back, PHTelem{numPatches}(i).nodesGlobal(back_nodes)];
        end
    end
end


%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);
bcdof_front = unique(bcdof_front);
bcdof_back = unique(bcdof_back);

%take into account that there are 3 global indices (displacements) per node
bcdof_all = [bcdof_right, bcdof_left, bcdof_front, bcdof_back];
bcdof = [3*bcdof_all-2, 3*bcdof_all-1, 3*bcdof_all];

stiff(bcdof,:)=[];  % zero out the rows and  columns of the K matrix
stiff(:,bcdof) = [];
mass(bcdof,:) = [];
mass(:, bcdof) = [];