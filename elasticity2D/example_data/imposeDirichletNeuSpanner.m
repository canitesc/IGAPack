function [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuSpanner(stiff, rhs, PHTelem, GIFTmesh, p, q, tauy)
%impose Dirichlet boundary conditions for spanner
%fixed (homogeneous) point boundary condtions at the corners
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

toleq = 1e-10;

%set the boundary degree of freedom and elements from the 1st patch
for indexPatch=1:6
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if indexPatch==6 && isempty(PHTelem{6}(i).neighbor_left) && abs(PHTelem{6}(i).vertex(2))<toleq
                bcdof_octop = PHTelem{6}(i).nodesGlobal(1);
            end
            if indexPatch==5 && isempty(PHTelem{5}(i).neighbor_left) && abs(PHTelem{5}(i).vertex(2))<toleq
                bcdof_ictop = PHTelem{5}(i).nodesGlobal(1);
            end
            if indexPatch==3 && isempty(PHTelem{3}(i).neighbor_left) && abs(PHTelem{3}(i).vertex(2))<toleq
                bcdof_icbot = PHTelem{3}(i).nodesGlobal(1);
            end
            if indexPatch==2 && isempty(PHTelem{2}(i).neighbor_left) && abs(PHTelem{2}(i).vertex(2))<toleq
                bcdof_ocbot = PHTelem{2}(i).nodesGlobal(1);
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_right)
                neumann_right = [neumann_right; i, 2, indexPatch];
            end
        end
    end
end


%impose the constraints (zero x and y displacement) in the linear system
bcdof = [2*bcdof_octop, 2*bcdof_ictop-1, 2*bcdof_ictop, 2*bcdof_icbot-1, 2*bcdof_icbot, 2*bcdof_ocbot];
bcval = zeros(size(bcdof));



%impose Neumann boundary conditons
neumann = neumann_right;
for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    patchIndex = neumann(i_neu,3);
    corient = neumann(i_neu, 2);
    
    xmin = PHTelem{patchIndex}(i).vertex(1);
    xmax = PHTelem{patchIndex}(i).vertex(3);
    ymin = PHTelem{patchIndex}(i).vertex(2);
    ymax = PHTelem{patchIndex}(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem{patchIndex}(i).C,1);
    
    
    scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
    
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
        R = phtBasis(u_hat, v_hat, PHTelem{patchIndex}(i).C, p, q);
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space
        [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, u_hat, v_hat, xmin, ymin, xmax, ymax);
        
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
        
        taux = 0;
        
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
