function [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletCubeWithHoleC0MP(stiff, rhs, PHTelem, GIFTmesh, p, q, r)
%impose Dirichlet and Neumann boundary conditions for the cube with a hole
%problem
%symmetry conditions on the boundaries that interesct the hole, traction in
%the vertical direction
% for C0 multipatch discretization


%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_front = [];
bcdof_right = [];
bcdof_down = [];
bcdof_back = [];
bcdof_up = [];

%define side node indices
down_nodes = 1:(p+1)*(q+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);
front_nodes = [];
neumann_down = [];
neumann_up = [];

numPatches = length(PHTelem);

for i=1:r+1
    front_nodes = [front_nodes, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end
%front_nodes
back_nodes = front_nodes + (p+1)*q;

for indexPatch=1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_front) && indexPatch<=2
                %x fixed in patches 1 and 2, front side                
                bcdof_front = [bcdof_front, PHTelem{indexPatch}(i).nodesGlobal(front_nodes)];
                
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_down)
                neumann_down = [neumann_down; i, 5, indexPatch];
                bcdof_down = [bcdof_down, PHTelem{indexPatch}(i).nodesGlobal(down_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                neumann_up = [neumann_up; i, 6, indexPatch];
                bcdof_up = [bcdof_up, PHTelem{indexPatch}(i).nodesGlobal(up_nodes)];
            end
            
            if isempty(PHTelem{indexPatch}(i).neighbor_left) && ((indexPatch == 1) || (indexPatch == 3))
                %z fixed in patches 1 and 3, left side
                bcdof_left = [bcdof_left, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_right) && ((indexPatch == 2) || (indexPatch == 4))
                %y fixed in patches 2 and 4, right side
                bcdof_right = [bcdof_right, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_back) && (indexPatch >= 3)
                %axis of revolution, fixed in the y and z directions in
                %patches 3 and 4
                bcdof_back = [bcdof_back, PHTelem{indexPatch}(i).nodesGlobal(back_nodes)];
            end
        end
    end
end

%remove duplicated entries
bcdof_front = unique(bcdof_front);
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);
bcdof_back = unique(bcdof_back);

%bcdof_down = unique(bcdof_down);

%take into account that there are 3 global indices (displacements) per node
bcdof_left_x = 3*bcdof_left-2;
bcdof_left_y = 3*bcdof_left-1;
bcdof_left_z = 3*bcdof_left;

bcdof_right_x = 3*bcdof_right-2;
bcdof_right_y = 3*bcdof_right-1;
bcdof_right_z = 3*bcdof_right;

bcdof_front_x = 3*bcdof_front-2;
bcdof_front_y = 3*bcdof_front-1;
bcdof_front_z = 3*bcdof_front;

bcdof_up_x = 3*bcdof_up-2;
bcdof_up_y = 3*bcdof_up-1;
bcdof_up_z = 3*bcdof_up;

bcdof_back_x = 3*bcdof_back-2;
bcdof_back_y = 3*bcdof_back-1;
bcdof_back_z = 3*bcdof_back;

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
%bcval_right_y = zeros(size(bcdof_right_y));
%bcval_right_z = zeros(size(bcdof_right_z));

% size(bcdof_left_y)
% size(bcdof_front_z)
% size(bcdof_right_x)
% size(bcdof_back_x)
% size(bcdof_back_y)
%
bcdof = unique([bcdof_left_z, bcdof_right_x, bcdof_front_y, bcdof_back_x, bcdof_back_z]);
%bcval = [bcval_left_y, bcval_front_z, bcval_right_x];
bcval = zeros(size(bcdof));
%impose Neumann boundary conditons
neumann = neumann_up;

%surface_area = 0;

ngauss_edge = p+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );
surface_area = 0;
for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    patchIndex = neumann(i_neu, 3);
    
    xmin = PHTelem{patchIndex}(i).vertex(1);
    xmax = PHTelem{patchIndex}(i).vertex(4);
    ymin = PHTelem{patchIndex}(i).vertex(2);
    ymax = PHTelem{patchIndex}(i).vertex(5);
    zmin = PHTelem{patchIndex}(i).vertex(3);
    zmax = PHTelem{patchIndex}(i).vertex(6);
    
    if (corient == 5) || (corient == 6)
        scalefac = (xmax-xmin)*(ymax-ymin)/4;
    elseif (corient == 2) || (corient == 4)
        scalefac = (ymax-ymin)*(zmax-zmin)/4;
    elseif (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)*(zmax-zmin)/4;
    end
    
    nument = size(PHTelem{patchIndex}(i).C,1);
    
    scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
    
    tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
    localrhsed = zeros(3*nument, 1);
    
    %loop over Gauss points and compute the integral
    for jgauss = 1:ngauss_edge
        for igauss =1:ngauss_edge
            
            %find the evaluation point, taking into account the side
            %orientation
            switch corient
                case 1
                    u_hat = gauss_coord_edge(igauss);
                    v_hat = -1;
                    w_hat = gauss_coord_edge(jgauss);
                case 2
                    u_hat = 1;
                    v_hat = gauss_coord_edge(igauss);
                    w_hat = gauss_coord_edge(jgauss);
                case 3
                    u_hat = gauss_coord_edge(igauss);
                    v_hat = 1;
                    w_hat = gauss_coord_edge(jgauss);
                case 4
                    u_hat = -1;
                    v_hat = gauss_coord_edge(igauss);
                    w_hat = gauss_coord_edge(jgauss);
                case 5
                    u_hat = gauss_coord_edge(igauss);
                    v_hat = gauss_coord_edge(jgauss);
                    w_hat = -1;
                case 6
                    u_hat = gauss_coord_edge(igauss);
                    v_hat = gauss_coord_edge(jgauss);
                    w_hat = 1;
            end
            
            %evaluate the basis functions
            R = phtBasis(u_hat, v_hat, w_hat, PHTelem{patchIndex}(i).C, p, q, r);
            
            %evaluate the derivatives of the mapping from parameter
            %space to physical space
            [coord, dxdxi] = paramMap3D( GIFTmesh{patchIndex}, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
            dxdxi = dxdxi';
             R = R';
%                                      plot3(coord(1),coord(2),coord(3),'r*')
%                                     drawnow
%                                     hold on
            %------------------------------------------------------------------------;
            %  computation of normal and jacobian;
            switch corient
                case 1
                    nor(1) = dxdxi(2,1)*dxdxi(3,3) - dxdxi(3,1)*dxdxi(2,3);
                    nor(2) = dxdxi(3,1)*dxdxi(1,3) - dxdxi(1,1)*dxdxi(3,3);
                    nor(3) = dxdxi(1,1)*dxdxi(2,3) - dxdxi(2,1)*dxdxi(1,3);
                case 2
                    nor(1) = dxdxi(2,2)*dxdxi(3,3) - dxdxi(3,2)*dxdxi(2,3);
                    nor(2) = dxdxi(3,2)*dxdxi(1,3) - dxdxi(1,2)*dxdxi(3,3);
                    nor(3) = dxdxi(1,2)*dxdxi(2,3) - dxdxi(2,2)*dxdxi(1,3);
                case 3
                    nor(1) = -dxdxi(2,1)*dxdxi(3,3) + dxdxi(3,1)*dxdxi(2,3);
                    nor(2) = -dxdxi(3,1)*dxdxi(1,3) + dxdxi(1,1)*dxdxi(3,3);
                    nor(3) = -dxdxi(1,1)*dxdxi(2,3) + dxdxi(2,1)*dxdxi(1,3);
                case 4
                    nor(1) = -dxdxi(2,2)*dxdxi(3,3) + dxdxi(3,2)*dxdxi(2,3);
                    nor(2) = -dxdxi(3,2)*dxdxi(1,3) + dxdxi(1,2)*dxdxi(3,3);
                    nor(3) = -dxdxi(1,2)*dxdxi(2,3) + dxdxi(2,2)*dxdxi(1,3);
                case 5
                    nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1);
                    nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1);
                    nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1);
                case 6
                    nor(1) = -dxdxi(2,2)*dxdxi(3,1) + dxdxi(3,2)*dxdxi(2,1);
                    nor(2) = -dxdxi(3,2)*dxdxi(1,1) + dxdxi(1,2)*dxdxi(3,1);
                    nor(3) = -dxdxi(1,2)*dxdxi(2,1) + dxdxi(2,2)*dxdxi(1,1);
            end
            
            J = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
            normal = nor/J; % normal vector in three dimensions
            
            [ ~, stress ] = exactSolSphericalCubeHole( coord );
            
            % R
            taux = normal(1)*stress(1) + normal(2)*stress(4)+ normal(3)*stress(6);
            tauy = normal(1)*stress(4) + normal(2)*stress(2)+ normal(3)*stress(5);
            tauz = normal(1)*stress(6) + normal(2)*stress(5)+ normal(3)*stress(3);
            
            localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(3:3:end) = localrhsed(3:3:end) + R.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            
            surface_area = surface_area + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
        end
    end
    rhs(tscrtx)=rhs(tscrtx)+localrhsed;
end

surface_area
[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);