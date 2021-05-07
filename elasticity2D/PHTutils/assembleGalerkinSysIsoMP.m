function [ stiff, rhs ] = assembleGalerkinSysIsoMP( PHTelem, controlPts, sizeBasis, p, q, Cmat )
%assembles the stiffness matrix and rhs (Galerkin method)
%uses GIFT mapping
%supports multipatches

%Gauss points
ngauss_x = p+1;
ngauss_y = q+1;

[gauss_weight_x, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[gauss_weight_y, gauss_coord_y] = quadrature( ngauss_y, 'GAUSS', 1 );

%take the transpose so that they are in the format expected by
%bernstein_basis
gauss_coord_x = gauss_coord_x';
gauss_coord_y = gauss_coord_y';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v] = bernstein_basis(gauss_coord_y,q);

dBdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
R = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));

%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
        R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
    end
end

%allocate memory for the triplet arrays
indexCounter = 0;
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nument = size(PHTelem{patchIndex}(i).C,1);
            indexCounter = indexCounter + 4*nument^2;
        end
    end
end

II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);
indexCounter = 0;

%initialize LHS stiffness matrix and RHS vector
dim = 2; %the dimension of physical space
rhs = zeros(dim*sizeBasis,1);

%assemble the stiffness matrix and RHS
elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter + 1;
            
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            nodes = PHTelem{indexPatch}(i).nodes(1:nument);
            dsctrx = reshape([2*sctrx-1; 2*sctrx],1,2*nument);
            
            localstiff = zeros(2*nument, 2*nument); %local stiffness
            cpts = controlPts{indexPatch}(nodes, 1:2);
            wgts = controlPts{indexPatch}(nodes, 3);
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    
                    dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                    dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                    RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,:));
                    RR = RR .* wgts;
                    dRdx = dRdx .* wgts;
                    dRdy = dRdy .* wgts;
                    
                    w_sum = sum(RR);
                    dw_xi = sum(dRdx);
                    dw_eta = sum(dRdy);
                    
                    dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                    dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    
                    dR  = [dRdx';dRdy'];
                    dxdxi = dR*cpts;
                    
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\dR;
                    J = det(dxdxi);
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dR(1,:);
                    B(2:2:2*nument,2) = dR(2,:);
                    B(1:2:2*nument-1,3) = dR(2,:);
                    B(2:2:2*nument,3) = dR(1,:);
                    
                    %TODO: implement non-zero volume force
                    localstiff = localstiff + B * Cmat * B' * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                end
            end
            II(indexCounter+1:indexCounter+4*nument^2) = repmat(dsctrx,1,2*nument);
            JJ(indexCounter+1:indexCounter+4*nument^2) = reshape(repmat(dsctrx,2*nument,1),1,4*nument^2);
            S(indexCounter+1:indexCounter+4*nument^2) = reshape(localstiff,1,4*nument^2);
            indexCounter = indexCounter + 4*nument^2;
        end
    end
end
stiff = sparse(II,JJ,S,2*sizeBasis,2*sizeBasis);
disp(['The mesh has ', num2str(elementCounter), ' active elements.'])

