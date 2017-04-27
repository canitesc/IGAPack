function [ stiff, rhs] = imposeDirichletNeuCantBeam(stiff, rhs, PHTelem, GIFTmesh, p, q, L, D, bound_press)
%impose Neuamann and Dirichlet boundary conditions for the L-shaped domain
%fix displacements in both directions at the re-entrant corner and in x
%direction at the outer corner
%supports multipatches
%produces a symmetric linear system

ngauss_edge = 2*q;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann =[];

neu_min = 13/16;
neu_max = 15/16;

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);




bcdof = [];

for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))                
                bcdof = [bcdof, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];                                
            end
            
            if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_up))                
                bcdof = [bcdof, PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];                                
            end
            
            if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_left))                
                bcdof = [bcdof, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];                                
            end
            
            if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_left))                
                bcdof = [bcdof, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];                                
            end
            
            if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_down)) 
                bcdof = [bcdof, PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];                                
            end
                        
            if (patchIndex==3) && isempty(PHTelem{patchIndex}(i).neighbor_up)
                %add elements that partially overlap on the left
                if (PHTelem{patchIndex}(i).vertex(1)<=neu_min) && (PHTelem{patchIndex}(i).vertex(3)>=neu_min)
                    neumann = [neumann; i, 3, patchIndex];
                end
                %add elements that are completely included in the traction
                %region
                if (PHTelem{patchIndex}(i).vertex(1)>=neu_min) && (PHTelem{patchIndex}(i).vertex(3)<=neu_max)
                    neumann = [neumann; i, 3, patchIndex];                    
                end
                %add elements that partially overlap on the right
                if (PHTelem{patchIndex}(i).vertex(1)<=neu_max) && (PHTelem{patchIndex}(i).vertex(3)>=neu_max)
                    neumann = [neumann; i, 3, patchIndex];
                end
            end
           
        end
    end
end


%impose the constraints (zero x and y displacement) in the linear system
bcdof = unique(bcdof);
bcdof = [2*bcdof-1, 2*bcdof];
bcval = zeros(size(bcdof));

neumann 
%impose Neumann boundary conditons
num_neumann_elem = size(neumann,1);

for i_neu=1:num_neumann_elem
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    patchIndex = neumann(i_neu, 3);
    
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
        taux = 0;
%         coord
%         (D + neu_min*L)
%         (D + neu_max*L)
        if (coord(1) >= (D + neu_min*L)) && (coord(1) <= (D + neu_max*L))
            tauy = -bound_press;
        else
            tauy = 0;
        end
     %   pause
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
