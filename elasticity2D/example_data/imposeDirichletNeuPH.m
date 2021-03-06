function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPH(stiff, rhs, PHTelem, GIFTmesh, p, q, rad, tx)
%impose Dirichlet boundary conditions for elastic rectangle
%fixed (homogeneous) boundary conditions on the left side
%Neumann (traction) boundary conditions on the right side

ngauss_edge = q+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );



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

for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        if isempty(PHTelem(i).neighbor_down)
            bcdof_down = [bcdof_down, PHTelem(i).nodesGlobal(down_nodes)];
            neumann_down = [neumann_down; i, 1];
        end
        if isempty(PHTelem(i).neighbor_right)
            bcdof_right = [bcdof_right, PHTelem(i).nodesGlobal(right_nodes)];
            neumann_right = [neumann_right; i, 2];
        end
        if isempty(PHTelem(i).neighbor_up)
            bcdof_up = [bcdof_up, PHTelem(i).nodesGlobal(up_nodes)];
            neumann_up = [neumann_up; i, 3];
        end
        if isempty(PHTelem(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem(i).nodesGlobal(left_nodes)];
            neumann_left = [neumann_left; i, 4];
        end
        
        
    end
end

%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);
bcdof_up = unique(bcdof_up);
bcdof_down = unique(bcdof_down);

bcdof_left_x = 2*bcdof_left-1;
bcdof_left_y = 2*bcdof_left;

bcdof_right_x = 2*bcdof_right-1;
bcdof_right_y = 2*bcdof_right;

bcdof_up_x = 2*bcdof_up-1;
bcdof_up_y = 2*bcdof_up;

bcdof_down_x = 2*bcdof_down-1;
bcdof_down_y = 2*bcdof_down;

%impose the symmetry boundary conditions on the top and bottom edge
bcval_right_x = zeros(size(bcdof_right_x));
bcval_left_y = zeros(size(bcdof_left_y));

bcdof = [bcdof_right_x, bcdof_left_y];
bcval = [bcval_right_x, bcval_left_y];


%impose Neumann boundary conditons
neumann = neumann_up;
for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    
    xmin = PHTelem(i).vertex(1);
    xmax = PHTelem(i).vertex(3);
    ymin = PHTelem(i).vertex(2);
    ymax = PHTelem(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem(i).C,1);
    
    
    scrtx = PHTelem(i).nodesGlobal(1:nument);
    
    if corient==1
        scrtx = scrtx(down_nodes);
        nument_edge = p+1;
    elseif corient==2
        scrtx = scrtx(right_nodes);
        nument_edge = q+1;
    elseif corient==3
        scrtx = scrtx(up_nodes);
        nument_edge = p+1;
    else
        scrtx = scrtx(left_nodes);
        nument_edge = q+1;
    end
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument_edge);
    localrhsed = zeros(2*nument_edge, 1);
    
    %loop over Gauss points and compute the integral
    for igauss = 1:ngauss_edge
        
        %find the evaluation point, taking into account the edge type
        if corient==1
            v_hat = -1;
            u_hat = gauss_coord_edge(igauss);
        elseif corient==2
            u_hat = 1;
            v_hat = gauss_coord_edge(igauss);
        elseif corient==3
            v_hat = 1;
            u_hat = gauss_coord_edge(igauss);
        else
            u_hat = -1;
            v_hat = gauss_coord_edge(igauss);
        end
        
        %evaluate the basis functions
        R = phtBasis(u_hat, v_hat, PHTelem(i).C, p, q);
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space
        [coord, dxdxi] = paramMap( GIFTmesh, u_hat, v_hat, xmin, ymin, xmax, ymax);        
        stress = ghole(coord(1), coord(2), rad, tx);
        dxdxi = dxdxi';
        %  jacobian of edge mapping;
        if((corient==1)||(corient==3))
            ee = dxdxi(1,1)^2+dxdxi(2,1)^2;
        else
            ee = dxdxi(1,2)^2+dxdxi(2,2)^2;
        end
        
        % Jacobian of face mapping
        J = sqrt(ee);%Length of edge
        %------------------------------------------------------------------------;
        %  computation of normal;
        if(corient==1)
            nor(1) = dxdxi(2,1);% dy/dxi
            nor(2) = -dxdxi(1,1);%dx/dxi
        elseif(corient==2)
            nor(1) = dxdxi(2,2);
            nor(2) = -dxdxi(1,2);
        elseif(corient==3)
            nor(1) = -dxdxi(2,1);
            nor(2) = dxdxi(1,1);
        else
            nor(1) = -dxdxi(2,2);
            nor(2) = dxdxi(1,2);
        end
        
        tmp = sqrt(nor(1)^2 + nor(2)^2);
        normal = nor/tmp; % normal vector in two dimensions
        
        
        %consider on the values corresponding to the shape functions on the
        %boundary
        if corient==1
            R = R(down_nodes)';
        elseif corient==2
            R = R(right_nodes)';
        elseif corient==3
            R = R(up_nodes)';
        else
            R = R(left_nodes)';
        end
        
        taux = normal(1)*stress(1) + normal(2)*stress(3);
        tauy = normal(1)*stress(3) + normal(2)*stress(2);
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
