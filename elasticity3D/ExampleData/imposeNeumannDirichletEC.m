function [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletEC(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_press)
%impose Dirichlet boundary conditions for elastic cube
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side

%detect which nodes have support on the left and right boundary
bcdof_top = [];
bcdof_right = [];

%define side node indices
down_nodes = 1:(p+1)*(q+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

neumann_down = [];

for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        if isempty(PHTelem(i).neighbor_down)
            %bcdof = [bcdof, PHTelem(i).nodes(down_nodes)];
            neumann_down = [neumann_down; i, 5];
        end
%         if isempty(PHTelem(i).neighbor_right)
%             bcdof_right = [bcdof_right, PHTelem(i).nodes(right_nodes)];
%         end
        if isempty(PHTelem(i).neighbor_up)
            bcdof_top = [bcdof_top, PHTelem(i).nodes(up_nodes)];
        end
%         if isempty(PHTelem(i).neighbor_left)
%             bcdof_left = [bcdof_left, PHTelem(i).nodes(left_nodes)];
%         end
    end
end

%remove duplicated entries
bcdof_top = unique(bcdof_top);
%bcdof_left = unique(bcdof_left);

%take into account that there are 3 global indices (displacements) per node
%bcdof_right_x = 3*bcdof_right-2;
%bcdof_right_y = 3*bcdof_right-1;
%bcdof_right_z = 3*bcdof_right;

bcdof_top_x = 3*bcdof_top-2;
bcdof_top_y = 3*bcdof_top-1;
bcdof_top_z = 3*bcdof_top;

%impose the prescribed displacement for each global index
bcval_top_x = zeros(size(bcdof_top_x));
bcval_top_y = zeros(size(bcdof_top_y));
bcval_top_z = zeros(size(bcdof_top_z));


%bcval_right_x = bound_disp*ones(size(bcdof_right_x));
%bcval_right_y = zeros(size(bcdof_right_y));
%bcval_right_z = zeros(size(bcdof_right_z));


bcdof = [bcdof_top_x, bcdof_top_y, bcdof_top_z];
bcval = [bcval_top_x, bcval_top_y, bcval_top_z];
%impose Neumann boundary conditons
neumann = neumann_down;
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
            
            %evaluate the basis functions
            R = phtBasis(u_hat, v_hat, w_hat, PHTelem(i).C, p, q, r);
            
            %evaluate the derivatives of the mapping from parameter
            %space to physical space
            [coord, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
            dxdxi = dxdxi';
%             plot3(coord(1),coord(2),coord(3),'r*')
%             drawnow
%             hold on
            %------------------------------------------------------------------------;
            %  computation of normal and jacobian;
            if(corient==5)
                nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1);
                nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1);
                nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1);
            end
            
            J = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
            normal = nor/J; % normal vector in two dimensions
            
            %consider on the values corresponding to the shape functions on the
            %boundary
            if corient==5
                R = R(down_nodes)';
            end
           
            
            %assume traction is constant
            taux = 0; %bound_press*normal(1)*ones(nument_side,1);
            tauy = -bound_press*ones(nument_side,1);
            tauz = 0; % bound_press*normal(3)*ones(nument_side,1);
            
            localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(3:3:end) = localrhsed(3:3:end) + R.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
        end
    end
   % localrhsed
    rhs(tscrtx)=rhs(tscrtx)+localrhsed;
   % pause
end
[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);