function [ rhs] = imposeNeumannIso(rhs, PHTelem, controlPts, neumann,...
                    traction_fun, p, q)
% Impose homoegeneous Dirichlet and Neumann boundary conditions 
% Inputs
% -------
%    rhs     - initial rhs 
%   PHTelem  - cell array of PHTelem structures containing the mesh information
% controlPts - cell array of control points for each patch
%  neumann   - array containing the patch indices and edge indices for the 
%                Neumann boundary conditions in the format [patchIndex,
%                 edgeIndex], where edgeIndex has the encoding 1-down(v=0),
%                 2 - right (u=1), 3-up (v=1),4-left(u=0)
% traction_fun - function handle of the type traction_fun(coords, normal, params)
%                 that returns the traction at points coord with the given 
%                 outer normal vector the the boundary, and params is a
%                 structure that contains additional information about the
%                 problem (e.g. material properties or geometry
%                 information)
%                on the Neumann boundary
%   p        - polynomial degree in the u direction
%   q        - polynomial degree in the v direction
%
%  Outputs
%  -------
%     rhs  - updated rhs vector


%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann_elem = [];
ngauss_edge = 2*max(p,q);
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

% loop over all patches and elements and check that they are on the Neumann
% or Dirichlet boundary
for indexPatch=1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            % check for the Dirichlet nodes
            for j = 1:size(neumann, 1)
                if indexPatch==neumann(j,1)
                    if isempty(PHTelem{indexPatch}(i).neighbor_down) && neumann(j,2)==1
                        neumann_elem = [neumann_elem; i, 1, indexPatch];
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_right) && neumann(j,2)==2
                        neumann_elem = [neumann_elem; i, 2, indexPatch];
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_up) && neumann(j,2)==3
                        neumann_elem = [neumann_elem; i, 3, indexPatch];
                    end
                    if isempty(PHTelem{indexPatch}(i).neighbor_left) && neumann(j,2)==4
                        neumann_elem = [neumann_elem; i, 4, indexPatch];
                    end
                end
            end
        end
    end
end


%impose Neumann boundary conditons
num_neumann_elem = size(neumann_elem,1);
boundaryLength = 0;
for i_neu=1:num_neumann_elem
    
    i = neumann_elem(i_neu,1);
    corient = neumann_elem(i_neu, 2);
    indexPatch = neumann_elem(i_neu, 3);
    
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
    
    sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
    
    if corient==1
        sctrx = sctrx(down_nodes);
        nument_edge = p+1;
    elseif corient==2
        sctrx = sctrx(right_nodes);
        nument_edge = q+1;
    elseif corient==3
        sctrx = sctrx(up_nodes);
        nument_edge = p+1;
    else
        sctrx = sctrx(left_nodes);
        nument_edge = q+1;
    end
    dsctrx = reshape([2*sctrx-1; 2*sctrx],1,2*nument_edge);
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
        hold on

        plot(coord(1), coord(2), '.r')
        
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
        
        hold on        
        plot(coord(1), coord(2),'.r')
                
        [taux, tauy] = traction_fun(coord, normal);
        boundaryLength = boundaryLength + scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dsctrx)=rhs(dsctrx)+localrhsed;
end
disp(['Neumann boundary length is ', num2str(boundaryLength)])


