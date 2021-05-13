function [ stiff, rhs] = imposeHomogeneousDirichlet(stiff, rhs, PHTelem,...
                        dirichlet, p, q)
% Impose homoegeneous Dirichlet boundary conditions 
% Inputs
% -------
%    stiff   - stiffness matrix
%    rhs     - rhs 
%   PHTelem  - cell array of PHTelem structures containing the mesh information
%  dirichlet - array containing the patch indices and edge indices for 
%              the homogeneous Dirichlet boundary conditions as a 2D matrix 
%              with each row of the form [patchIndex, edgeIndex, fix_x, fix_y]
%              where the edge index uses the usual encoding: 1 - down (v=0), 
%              2 - right (u=1),  3 - up (v=1), 4 - left (u=0) and the 
%              columns fix_x and fix_y are either 0 or 1 depending on
%              whether the x or the y direction (or both) are fixed
%   p        - polynomial degree in the u direction
%   q        - polynomial degree in the v direction
%
%  Outputs
%  -------
%    stiff - updated stiffness matrix
%     rhs  - updated rhs vector

bcdof_x = [];
bcdof_y = [];


%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

% loop over all patches and elements and check that they are on the Neumann
% or Dirichlet boundary
for indexPatch=1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            % check for the Dirichlet nodes
            for j = 1:size(dirichlet, 1)
                if indexPatch==dirichlet(j,1)
                    if isempty(PHTelem{indexPatch}(i).neighbor_down) && dirichlet(j,2)==1
                        if dirichlet(j,3)==1
                            bcdof_x = [bcdof_x, PHTelem{indexPatch}(i).nodesGlobal(down_nodes)];
                        end
                        if dirichlet(j,4)==1
                            bcdof_y = [bcdof_y, PHTelem{indexPatch}(i).nodesGlobal(down_nodes)];
                        end
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_right) && dirichlet(j,2)==2
                        if dirichlet(j,3)==1
                            bcdof_x = [bcdof_x, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
                        end
                        if dirichlet(j,4)==1
                            bcdof_y = [bcdof_y, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
                        end
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_up) && dirichlet(j,2)==3
                        if dirichlet(j,3)==1
                            bcdof_x = [bcdof_x, PHTelem{indexPatch}(i).nodesGlobal(up_nodes)];
                        end
                        if dirichlet(j,4)==1
                            bcdof_y = [bcdof_y, PHTelem{indexPatch}(i).nodesGlobal(up_nodes)];
                        end
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_left) && dirichlet(j,2)==4
                        if dirichlet(j,3)==1
                            bcdof_x = [bcdof_x, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
                        end
                        if dirichlet(j,4)==1
                            bcdof_y = [bcdof_y, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
                        end
                    end
                end
            end
        end
    end
end

bcdof = unique([2*bcdof_x-1, 2*bcdof_y]);
bcval = zeros(size(bcdof));

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
