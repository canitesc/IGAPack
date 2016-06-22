function [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletCubeWithHole(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_trac)
%impose Dirichlet and Neumann boundary conditions for the cube with a hole
%problem
%symmetry conditions on the boundaries that interesct the hole, traction in
%the vertical direction


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

for i=1:r+1
    front_nodes = [front_nodes, (p+1)*(q+1)*(i-1)+1:(p+1)*((q+1)*(i-1)+1)];
end
%front_nodes
back_nodes = front_nodes + (p+1)*q;


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
            neumann_up = [neumann_up; i, 6];
            bcdof_up = [bcdof_up, PHTelem(i).nodes(up_nodes)];
        end
        
        if isempty(PHTelem(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem(i).nodes(left_nodes)];
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
bcval_left_x = zeros(size(bcdof_left_x));
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

for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    
    xmin = PHTelem(i).vertex(1);
    xmax = PHTelem(i).vertex(4);
    ymin = PHTelem(i).vertex(2);
    ymax = PHTelem(i).vertex(5);
    zmin = PHTelem(i).vertex(3);
    zmax = PHTelem(i).vertex(6);
    
    if (corient == 5) || (corient == 6)
        scalefac = (xmax-xmin)*(ymax-ymin)/4;
    end
    
    nument = size(PHTelem(i).C,1);
    
    scrtx = PHTelem(i).nodes(1:nument);
    
    if corient==5
        scrtx = scrtx(down_nodes);
        nument_side = (p+1)*(q+1);
    end
    
    if corient==6
        scrtx = scrtx(up_nodes);
        nument_side = (p+1)*(q+1);
    end
    
    tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument_side);
    localrhsed = zeros(3*nument_side, 1);
    
    %loop over Gauss points and compute the integral
    for jgauss = 1:ngauss_edge
        for igauss =1:ngauss_edge
            
            %find the evaluation point, taking into account the side
            %orientation
            if corient==5
                u_hat = gauss_coord_edge(igauss);
                v_hat = gauss_coord_edge(jgauss);
                w_hat = -1;
            end
            
             if corient==6
                u_hat = gauss_coord_edge(igauss);
                v_hat = gauss_coord_edge(jgauss);
                w_hat = 1;
            end
                       
            %evaluate the basis functions
            R = phtBasis(u_hat, v_hat, w_hat, PHTelem(i).C, p, q, r);
            
            %evaluate the derivatives of the mapping from parameter
            %space to physical space
            [coord, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
            dxdxi = dxdxi';
            
            [ disp, stress ] = exactSolSphericalCubeHole( coord );
%              plot3(coord(1),coord(2),coord(3),'r*')
%             drawnow
%             hold on
            %------------------------------------------------------------------------;
            %  computation of normal and jacobian;
            if(corient==5)
               nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1);
               nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1);
               nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1);     
            end
            
             if(corient==6)
               nor(1) = -dxdxi(2,2)*dxdxi(3,1) + dxdxi(3,2)*dxdxi(2,1);
               nor(2) = -dxdxi(3,2)*dxdxi(1,1) + dxdxi(1,2)*dxdxi(3,1);
               nor(3) = -dxdxi(1,2)*dxdxi(2,1) + dxdxi(2,2)*dxdxi(1,1);     
            end
            
            J = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
            normal = nor/J; % normal vector in three dimensions
                
            %consider on the values corresponding to the shape functions on the
            %boundary
            if corient==5
                R = R(down_nodes)';    
            end
            
             if corient==6
                R = R(up_nodes)';    
            end
            
           % R
            %assume traction is constant
            
           % plot3(coord(1),coord(2),coord(3),'.r')
           % hold on
            taux = normal(1)*stress(1) + normal(2)*stress(4)+ normal(3)*stress(6);
            tauy = normal(1)*stress(4) + normal(2)*stress(2)+ normal(3)*stress(5);
            tauz = normal(1)*stress(6) + normal(2)*stress(5)+ normal(3)*stress(3);
%            pause
%             taux = stress(1)*normal(1)*ones(nument_side,1);
%             tauy = stress(2)*normal(2)*ones(nument_side,1);
%             tauz = stress(3)*normal(3)*ones(nument_side,1);
            
            localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(3:3:end) = localrhsed(3:3:end) + R.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;   
            
            %surface_area = surface_area + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
        end
    end
    rhs(tscrtx)=rhs(tscrtx)+localrhsed;
end

%check that the surface area is computed correctly
%surface_area - pi/2

%rhs(sort([3*bcdof_down-2, 3*bcdof_down-1, 3*bcdof_down]))
%pause
[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);
