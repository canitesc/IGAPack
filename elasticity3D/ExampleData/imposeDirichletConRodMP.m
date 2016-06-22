function [ stiff, rhs, bcdof, bcval ] = imposeDirichletConRodMP(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_press)
%impose Dirichlet and Neumann boundary conditions for Connecting Rod 
%fixed (homogeneous) boundary conditions on the left side
%traction to the right on the right side

%detect which nodes have support on the left and right boundary
bcdof_fixed = [];
neumann_front =[];


numPatches = length(PHTelem);

%define side node indices
%down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
%up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

for indexPatch = 1:3
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_left)
                bcdof_fixed = [bcdof_fixed, PHTelem{indexPatch}(i).nodesGlobal(left_nodes)];
            end
        end
    end
end

for indexPatch = 7:9   
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_right)
                bcdof_fixed = [bcdof_fixed, PHTelem{indexPatch}(i).nodesGlobal(right_nodes)];
            end
        end
    end
end

%Neumann boundary conditions
for indexPatch = 11:25
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if isempty(PHTelem{indexPatch}(i).neighbor_front)
                neumann_front = [neumann_front; i, 1, indexPatch];
            end
        end
    end
end

%remove duplicated entries
bcdof_fixed = unique(bcdof_fixed);

%take into account that there are 3 global indices (displacements) per node
bcdof_fixed_x = 3*bcdof_fixed-2;
bcdof_fixed_y = 3*bcdof_fixed-1;
bcdof_fixed_z = 3*bcdof_fixed;


%impose the prescribed displacement for each global index
bcval_fixed_x = zeros(size(bcdof_fixed_x));
bcval_fixed_y = zeros(size(bcdof_fixed_y));
bcval_fixed_z = zeros(size(bcdof_fixed_z));


bcdof = [bcdof_fixed_x, bcdof_fixed_y, bcdof_fixed_z];
bcval = [bcval_fixed_x, bcval_fixed_y, bcval_fixed_z];


ngauss_edge = p+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );
neumann = neumann_front;
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
%                          plot3(coord(1),coord(2),coord(3),'r*')
%                         drawnow
%                         hold on
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
            %[patchIndex, i, J]
            normal = nor/J; % normal vector in three dimensions
            %assume traction is constant
            %stress = [0,0,-bound_press,0,0,0];
            
            %taux = normal(1)*stress(1) + normal(2)*stress(4)+ normal(3)*stress(6);
            %tauy = normal(1)*stress(4) + normal(2)*stress(2)+ normal(3)*stress(5);
            %tauz = normal(1)*stress(6) + normal(2)*stress(5)+ normal(3)*stress(3);
            taux = 0;
            tauy = 0;
            tauz = -bound_press;
            
            localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            localrhsed(3:3:end) = localrhsed(3:3:end) + R.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
            
            surface_area = surface_area + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
        end
    end
    rhs(tscrtx)=rhs(tscrtx)+localrhsed;
    %pause
end

surface_area

[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);