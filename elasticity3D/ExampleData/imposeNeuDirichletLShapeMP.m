function [ stiff, rhs, bcdof, bcval ] = imposeNeuDirichletLShapeMP(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_disp)
%impose Dirichlet boundary conditions for elastic crack
%fixed (homogeneous) boundary conditions on the bottom side
%displacement (bound_disp units) to the top on the top side

%detect which nodes have support on the left and right boundary
bcdof_bottom = [];
bcdof_top = [];

numPatches = length(PHTelem);

%for each neumann edge store the element index and orientation
%orientation: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up
neumann_up =[];
neumann_down =[];
neumann_left =[];
neumann_right =[];

bcdof_UpRightReentrant = [];
bcdof_UpRightOuter = [];

%define side node indices
down_nodes = 1:(p+1)*(q+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);
up_right_nodes = intersect(right_nodes, up_nodes);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);

for i=1:r+1
    front_nodes((i-1)*(p+1)+1:i*(p+1)) = (p+1)*(q+1)*(i-1)+1:(p+1)*(q+1)*(i-1)+p+1;
    back_nodes((i-1)*(p+1)+1:i*(p+1)) = i*(p+1)*(q+1):-1:i*(p+1)*(q+1)-p;
end

surface_area = 0;

% for patchIndex = 1:length(PHTelem)
%     for i=1:length(PHTelem{patchIndex})
%         if isempty(PHTelem{patchIndex}(i).children)
%             if (patchIndex==3) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
%                 neumann_up = [neumann_up; i, 6, patchIndex];
%             end
%             
%             if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
%                 neumann_up = [neumann_up; i, 6, patchIndex];
%             end
%             
%             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
%                 neumann_left = [neumann_left; i, 4, patchIndex];
%             end
%             
%             if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
%                 neumann_left = [neumann_left; i, 4, patchIndex];
%             end
%             
%             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
%                 neumann_right = [neumann_right; i, 2, patchIndex];
%             end
%             
%             if (patchIndex==3) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
%                 neumann_right = [neumann_right; i, 2, patchIndex];
%             end
%             
%             
%             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
%                 neumann_down = [neumann_down; i, 5, patchIndex];
%             end
%             
%             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
%                 neumann_down = [neumann_down; i, 5, patchIndex];
%             end
%             
%             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
%                 if (PHTelem{patchIndex}(i).vertex(4)==1) && (PHTelem{patchIndex}(i).vertex(6)==1)
%                     bcdof_UpRightReentrant = [bcdof_UpRightReentrant, PHTelem{patchIndex}(i).nodesGlobal(up_right_nodes)];
%                 end
%             end
%             
%             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
%                 if (PHTelem{patchIndex}(i).vertex(4)==1) && (PHTelem{patchIndex}(i).vertex(6)==1)
%                     bcdof_UpRightOuter = [bcdof_UpRightOuter, PHTelem{patchIndex}(i).nodesGlobal(up_right_nodes)];
%                 end
%             end
%             
%         end
%     end
% end
bcdof_left = [];
for i=1:length(PHTelem{16})
    if isempty(PHTelem{16}(i).children)        
        if isempty(PHTelem{16}(i).neighbor_back) && PHTelem{16}(i).vertex(5)==1
            bcdof_left = [bcdof_left, PHTelem{16}(i).nodesGlobal(back_nodes)];  
        end
     
    end
end

bcdof_right = [];
for i=1:length(PHTelem{4})
    if isempty(PHTelem{4}(i).children)        
        if isempty(PHTelem{4}(i).neighbor_right) && PHTelem{4}(i).vertex(4)==1
            bcdof_right = [bcdof_right, PHTelem{4}(i).nodesGlobal(right_nodes)];  
        end     
    end
end

%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

%take into account that there are 2 global indices (displacements) per node
bcdof_right_x = 3*bcdof_right-2;
bcdof_right_y = 3*bcdof_right-1;
bcdof_right_z = 3*bcdof_right;

bcdof_left_x = 3*bcdof_left-2;
bcdof_left_y = 3*bcdof_left-1;
bcdof_left_z = 3*bcdof_left;

%impose the prescribed displacement for each global index
bcval_left_x = zeros(size(bcdof_left_x));
bcval_left_y = zeros(size(bcdof_left_y));
bcval_left_z = zeros(size(bcdof_left_z));

bcval_right_x = zeros(size(bcdof_right_x));
bcval_right_y = -bound_disp*ones(size(bcdof_right_y));
bcval_right_z = bound_disp*ones(size(bcdof_right_z));

bcdof = [bcdof_left_x, bcdof_left_y, bcdof_left_z, bcdof_right_y];
bcval = [bcval_left_x, bcval_left_y, bcval_left_z, bcval_right_y];

% ngauss_edge = p+2;
% [gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );
% neumann = [neumann_up; neumann_down; neumann_left; neumann_right];


% for i_neu=1:size(neumann,1)
%     
%     i = neumann(i_neu,1);
%     corient = neumann(i_neu, 2);
%     patchIndex = neumann(i_neu, 3);
% 
%     xmin = PHTelem{patchIndex}(i).vertex(1);
%     xmax = PHTelem{patchIndex}(i).vertex(4);
%     ymin = PHTelem{patchIndex}(i).vertex(2);
%     ymax = PHTelem{patchIndex}(i).vertex(5);
%     zmin = PHTelem{patchIndex}(i).vertex(3);
%     zmax = PHTelem{patchIndex}(i).vertex(6);
%     
%     if (corient == 5) || (corient == 6)
%         scalefac = (xmax-xmin)*(ymax-ymin)/4;
%     elseif (corient == 2) || (corient == 4)
%         scalefac = (ymax-ymin)*(zmax-zmin)/4;
%     elseif (corient == 1) || (corient == 3)
%         scalefac = (xmax-xmin)*(zmax-zmin)/4;
%     end
%     
%     nument = size(PHTelem{patchIndex}(i).C,1);
%     
%     scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
%     
%     tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
%     localrhsed = zeros(3*nument, 1);
%     
%     %loop over Gauss points and compute the integral
%     for jgauss = 1:ngauss_edge
%         for igauss =1:ngauss_edge
%             
%             %find the evaluation point, taking into account the side
%             %orientation
%             switch corient
%                 case 2
%                     u_hat = 1;
%                     v_hat = gauss_coord_edge(igauss);
%                     w_hat = gauss_coord_edge(jgauss);
%                 case 4
%                     u_hat = -1;
%                     v_hat = gauss_coord_edge(igauss);
%                     w_hat = gauss_coord_edge(jgauss);
%                 case 5
%                     u_hat = gauss_coord_edge(igauss);
%                     v_hat = gauss_coord_edge(jgauss);
%                     w_hat = -1;
%                 case 6
%                     u_hat = gauss_coord_edge(igauss);
%                     v_hat = gauss_coord_edge(jgauss);
%                     w_hat = 1;
%             end
%             
%             %evaluate the basis functions          
%             R = phtBasis(u_hat, v_hat, w_hat, PHTelem{patchIndex}(i).C, p, q, r);
%             
%             %evaluate the derivatives of the mapping from parameter
%             %space to physical space
%             [coord, dxdxi] = paramMap3D( GIFTmesh{patchIndex}, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
%             dxdxi = dxdxi';
%             R = R';
% %                          plot3(coord(1),coord(2),coord(3),'r*')
% %                         drawnow
% %                         hold on
%             %------------------------------------------------------------------------;
%             %  computation of normal and jacobian;
%             switch corient
%                 case 1
%                     nor(1) = dxdxi(2,1)*dxdxi(3,3) - dxdxi(3,1)*dxdxi(2,3);
%                     nor(2) = dxdxi(3,1)*dxdxi(1,3) - dxdxi(1,1)*dxdxi(3,3);
%                     nor(3) = dxdxi(1,1)*dxdxi(2,3) - dxdxi(2,1)*dxdxi(1,3);
%                 case 2
%                     nor(1) = dxdxi(2,2)*dxdxi(3,3) - dxdxi(3,2)*dxdxi(2,3);
%                     nor(2) = dxdxi(3,2)*dxdxi(1,3) - dxdxi(1,2)*dxdxi(3,3);
%                     nor(3) = dxdxi(1,2)*dxdxi(2,3) - dxdxi(2,2)*dxdxi(1,3);
%                 case 3
%                     nor(1) = -dxdxi(2,1)*dxdxi(3,3) + dxdxi(3,1)*dxdxi(2,3);
%                     nor(2) = -dxdxi(3,1)*dxdxi(1,3) + dxdxi(1,1)*dxdxi(3,3);
%                     nor(3) = -dxdxi(1,1)*dxdxi(2,3) + dxdxi(2,1)*dxdxi(1,3);
%                 case 4
%                     nor(1) = -dxdxi(2,2)*dxdxi(3,3) + dxdxi(3,2)*dxdxi(2,3);
%                     nor(2) = -dxdxi(3,2)*dxdxi(1,3) + dxdxi(1,2)*dxdxi(3,3);
%                     nor(3) = -dxdxi(1,2)*dxdxi(2,3) + dxdxi(2,2)*dxdxi(1,3);
%                 case 5
%                     nor(1) = dxdxi(2,2)*dxdxi(3,1) - dxdxi(3,2)*dxdxi(2,1);
%                     nor(2) = dxdxi(3,2)*dxdxi(1,1) - dxdxi(1,2)*dxdxi(3,1);
%                     nor(3) = dxdxi(1,2)*dxdxi(2,1) - dxdxi(2,2)*dxdxi(1,1);
%                 case 6
%                     nor(1) = -dxdxi(2,2)*dxdxi(3,1) + dxdxi(3,2)*dxdxi(2,1);
%                     nor(2) = -dxdxi(3,2)*dxdxi(1,1) + dxdxi(1,2)*dxdxi(3,1);
%                     nor(3) = -dxdxi(1,2)*dxdxi(2,1) + dxdxi(2,2)*dxdxi(1,1);
%             end
%             
%             J = sqrt(nor(1)^2 + nor(2)^2 + nor(3)^2);
%             normal = nor/J; % normal vector in three dimensions
%             
%             [theta, rad] = cart2pol(coord(1),coord(3));
%             cracklength = 0.5;
%             Kexact1=0;
%             Kexact2 = KII;
%             stress=exact_stresses2(rad, theta, Kexact1, Kexact2);
%             
%             % R
%             taux = normal(1)*stress(1) + normal(2)*stress(4)+ normal(3)*stress(6);
%             tauy = normal(1)*stress(4) + normal(2)*stress(2)+ normal(3)*stress(5);
%             tauz = normal(1)*stress(6) + normal(2)*stress(5)+ normal(3)*stress(3);                        
%                       
%             localrhsed(1:3:end-2) = localrhsed(1:3:end-2) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
%             localrhsed(2:3:end-1) = localrhsed(2:3:end-1) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
%             localrhsed(3:3:end) = localrhsed(3:3:end) + R.*tauz.*scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
%             
%             surface_area = surface_area + scalefac.*gauss_weight_edge(igauss).*gauss_weight_edge(jgauss).*J;
%         end
%     end
%     rhs(tscrtx)=rhs(tscrtx)+localrhsed;
% end
% 
% surface_area
[stiff,rhs] = imposeDirichlet3d(stiff,rhs,bcdof,bcval);