function [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletCylinder(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_press)
%impose Dirichlet and Neumann boundary conditions for the hollow sphere
%symmetry condition on the flat edges, pressure on the interior surface

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
neumann_left = [];



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
            bcdof_up = [bcdof_up, PHTelem(i).nodes(up_nodes)];
        end
        
        if isempty(PHTelem(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem(i).nodes(left_nodes)];
            neumann_left = [neumann_left; i, 4];
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

bcdof_down = unique(bcdof_down);

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
bcdof = unique([bcdof_front_y, bcdof_back_x,bcdof_down_z,bcdof_up_z]);
%bcval = [bcval_left_y, bcval_front_z, bcval_right_x];
bcval = zeros(size(bcdof));

%impose Neumann boundary conditons
neumann = neumann_left;

ngauss_edge = p+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );
size(neumann,1)
norsurf = 0;
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
    if (corient == 4) || (corient == 2)
        scalefac = (ymax-ymin)*(zmax-zmin)/4;
    end
    
    nument = size(PHTelem(i).C,1);
    
    scrtx = PHTelem(i).nodes(1:nument)
    
    %     if corient==5
    %         scrtx = scrtx(down_nodes);
    %         nument_side = (p+1)*(q+1);
    %     elseif corient==4
    %         scrtx = scrtx(left_nodes);
    %         nument_side = (q+1)*(r+1);
    %     end
    tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
    %tscrtx = zeros(1,3*nument);
    localrhsed = zeros(3*nument, 1);
    
    %loop over Gauss points and compute the integral
    for jgauss = 1:ngauss_edge
        for igauss =1:ngauss_edge
            
            %find the evaluation point, taking into account the side
            %orientation
            if corient==5
                u_hat = gauss_coord_edge(igauss);
                v_hat = gauss_coord_edge(jgauss);
                w_hat = -1;
            elseif corient==4
                u_hat = -1;
                v_hat = gauss_coord_edge(igauss);
                w_hat = gauss_coord_edge(jgauss);
            end
          
            
            %evaluate the basis functions
            R = phtBasis(u_hat, v_hat, w_hat, PHTelem(i).C, p, q, r);
            
            
            %evaluate the derivatives of the mapping from parameter
            %space to physical space
            [coord, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
            coord
            dxdxi = dxdxi';
            %              plot3(coord(1),coord(2),coord(3),'r*')
            %             drawnow
            %             hold on
            %------------------------------------------------------------------------;
            %  computation of normal and jacobian;
            if(corient==5)
                nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1);
                nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1);
                nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1);
            elseif corient==4
                nor(1) = -dxdxi(2,2)*dxdxi(3,3) + dxdxi(3,2)*dxdxi(2,3);
                nor(2) = -dxdxi(3,2)*dxdxi(1,3) + dxdxi(1,2)*dxdxi(3,3);
                nor(3) = -dxdxi(1,2)*dxdxi(2,3) + dxdxi(2,2)*dxdxi(1,3);
            end
            
            
            J = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
            %coord
            %nor
            normal = nor/J; % normal vector in three dimensions
            %    pause
            %             %consider on the values corresponding to the shape functions on the
            %             %boundary
            %             if corient==5
            %                 R = R(down_nodes)';
            %             elseif corient==4
            %                 R = R(left_nodes)';
            %             end
            
            %% R
            % R(left_nodes)
            % pause

             norsurf = norsurf + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
             
            taux = -bound_press*normal(1)*ones(nument,1);
            tauy = -bound_press*normal(2)*ones(nument,1);
            tauz = -bound_press*normal(3)*ones(nument,1);         
            localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R'.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R'.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(3:3:end) = localrhsed(3:3:end) + R'.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
%             norsurf = norsurf + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
        end
        
    end
    
    rhs(tscrtx)=rhs(tscrtx)+localrhsed;
end
norsurf-pi/2
%rhs(sort([3*bcdof_down-2, 3*bcdof_down-1, 3*bcdof_down]))
%pause
rhs_old = rhs;
[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);
[rhs_old, rhs]
pause
