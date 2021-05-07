function [ stiff, rhs] = imposeDirichletNeuLSExSymIso(stiff, rhs, PHTelem, controlPts, p, q)
%impose Neuamann and Dirichlet boundary conditions for the L-shaped domain
%fix displacements in both directions at the re-entrant corner and in x
%direction at the outer corner
%supports multipatches
%produces a symmetric linear system

ngauss_edge = 2*q;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann_right = [];
neumann_up =[];


%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

for indexPatch=1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            if (indexPatch==2) && (isempty(PHTelem{indexPatch}(i).neighbor_left))
                if (PHTelem{indexPatch}(i).vertex(1)==0) && (PHTelem{indexPatch}(i).vertex(2)==0)
                    bcdof_reentrant = PHTelem{indexPatch}(i).nodesGlobal(1);                    
                end
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_right)
                neumann_right = [neumann_right; i, 2, indexPatch];
            end
            if isempty(PHTelem{indexPatch}(i).neighbor_up)
                neumann_up = [neumann_up; i, 3, indexPatch];
                if (indexPatch==2) && (PHTelem{indexPatch}(i).vertex(3)==1) && (PHTelem{indexPatch}(i).vertex(4)==1)
                    bcdof_outer = PHTelem{indexPatch}(i).nodesGlobal(end);                    
                end
            end            
        end
    end
end


%impose the constraints (zero x and y displacement) in the linear system
bcdof = [2*bcdof_reentrant-1, 2*bcdof_reentrant, 2*bcdof_outer];
bcval = zeros(size(bcdof));

%impose Neumann boundary conditons
neumann = [neumann_up; neumann_right];
num_neumann_elem = size(neumann,1);

for i_neu=1:num_neumann_elem
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    indexPatch = neumann(i_neu, 3);
    
    xmin = PHTelem{indexPatch}(i).vertex(1);
    xmax = PHTelem{indexPatch}(i).vertex(3);
    ymin = PHTelem{indexPatch}(i).vertex(2);
    ymax = PHTelem{indexPatch}(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem{indexPatch}(i).C,1);
    nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            
    cpts = controlPts{indexPatch}(nodes, 1:2);
    wgts = controlPts{indexPatch}(nodes, 3);
    
    
    scrtx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
    
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
        [R, dR] = phtBasisIso(u_hat, v_hat, PHTelem{indexPatch}(i).C, p, q, wgts);
        dR(1,:) = dR(1,:)*2/(xmax-xmin);
        dR(2,:) = dR(2,:)*2/(ymax-ymin);
        w_sum = sum(R);
        dw_xi = sum(dR(1,:));
        dw_eta = sum(dR(2,:));
                    
        dR(1,:) = dR(1,:)/w_sum - R*dw_xi/w_sum^2;
        dR(2,:) = dR(2,:)/w_sum - R*dw_eta/w_sum^2;
        R = R/w_sum;
                
        coord = R*cpts;
        
        dxdxi = dR*cpts;
        dxdxi = dxdxi';
        %  jacobian of edge mapping;
        if((corient==1)||(corient==3))
            J = hypot(dxdxi(1,1), dxdxi(2,1));
        else
            J = hypot(dxdxi(1,2), dxdxi(2,2));
        end
        
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
        
        tmp = hypot(nor(1),nor(2));
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
        
        [theta, r] = cart2pol(coord(1),coord(2));
        
        exactStress=getWedgeStresses( r, theta );
        
        taux = normal(1)*exactStress(1) + normal(2)*exactStress(3);
        tauy = normal(1)*exactStress(3) + normal(2)*exactStress(2);
                
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
