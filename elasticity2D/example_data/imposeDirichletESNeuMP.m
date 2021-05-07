function [ stiff, rhs, bcdof, bcval ] = imposeDirichletESNeuMP(stiff, rhs, PHTelem, GIFTmesh, p, q, bound_trac)
%impose Dirichlet boundary conditions for elastic rectangle
%fixed (homogeneous) boundary conditions on the left side
%Neumann (traction) boundary conditions on the right side
%supports multipatches

ngauss_edge = q+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

numPatches = length(PHTelem);


%detect which nodes have support on the left and right boundary
bcdof_left = [];
%bcdof_right = [];

%for each neumann edge store the element index and orientation
%orientation: 1-bottom, 2-right, 3-top, 4-left
neumann_right = [];

%define side node indices
%down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
%up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

%find the nodes corresponding to the Dirichlet boundary on patch 1
for i=1:length(PHTelem{1})
    if isempty(PHTelem{1}(i).children)
        if isempty(PHTelem{1}(i).neighbor_left)
            bcdof_left = [bcdof_left, PHTelem{1}(i).nodesGlobal(left_nodes)];
        end                
    end
end

%find the elements corresponding to the Neumann boundary on the last patch
for i=1:length(PHTelem{numPatches})
    if isempty(PHTelem{numPatches}(i).children)        
        if isempty(PHTelem{numPatches}(i).neighbor_right)            
            neumann_right = [neumann_right; i, 2];
        end    
    end
end

%remove duplicated entries
%bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);

bcdof_left_x = 2*bcdof_left-1;
bcdof_left_y = 2*bcdof_left;

%impose the prescribed displacement for each global index
bcval_left_x = zeros(size(bcdof_left_x));
bcval_left_y = zeros(size(bcdof_left_y));

bcdof = [bcdof_left_x, bcdof_left_y];
bcval = [bcval_left_x, bcval_left_y];


%impose Neumann boundary conditons
neumann = neumann_right;
for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    
    xmin = PHTelem{numPatches}(i).vertex(1);
    xmax = PHTelem{numPatches}(i).vertex(3);
    ymin = PHTelem{numPatches}(i).vertex(2);
    ymax = PHTelem{numPatches}(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem{numPatches}(i).C,1);
    nument_edge = q+1;
    localrhsed = zeros(2*nument_edge, 1);
    scrtx = PHTelem{numPatches}(i).nodesGlobal(1:nument);
    %assume corient = 2!
    scrtx = scrtx(right_nodes);
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument_edge);
    
    
    %loop over Gauss points and compute the integral
    for igauss = 1:ngauss_edge
        
        %evaluate the basis functions, assume corient = 2;
        u_hat = 1;
        v_hat = gauss_coord_edge(igauss);
        
        R = phtBasis(u_hat, v_hat, PHTelem{numPatches}(i).C, p, q);
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space, assume corient=2
        [~, dxdxi] = paramMap( GIFTmesh{numPatches}, u_hat, v_hat, xmin, ymin, xmax, ymax);
     %   dxdxi = dxdxi';
        J = sqrt(dxdxi(2,1)^2+dxdxi(2,2)^2);
        
        %consider on the values corresponding to the shape functions on the
        %boundary
        R = R(right_nodes)';
        
        %assume traction is constant
        taux = bound_trac(1)*ones(nument_edge,1);
        tauy = bound_trac(2)*ones(nument_edge,1);
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2(stiff, rhs, bcdof, bcval);
