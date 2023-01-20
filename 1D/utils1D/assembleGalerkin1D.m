function [ stiff, rhs, elementCounter ] = assembleGalerkin1D( PHTelem, GeoMesh, sizeBasis, p )
%assembles the stiffness matrix and rhs (Galerkin method)
%uses GeoMesh mapping
%supports multipatches

%Gauss points
if p>3
    ngauss_x = p+1;
else
    ngauss_x = p+2;
end
[gaussWeightU, gaussCoordU] = quadrature( ngauss_x, 'GAUSS', 1 );

%take the transpose so that they are in the format expected by
%bernstein_basis
gaussCoordU = gaussCoordU';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[Bu, dBdu] = bernstein_basis(gaussCoordU,p);

%allocate memory for the triplet arrays
indexCounter = 0;
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nument = size(PHTelem{patchIndex}(i).C,1);
            indexCounter = indexCounter + nument^2;
        end
    end
end
%indexCounter
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);

indexCounter = 0;

%initialize LHS stiffness matrix and RHS vector
dim = 1; %the dimension of physical space
%stiff = sparse(dim*sizeBasis,dim*sizeBasis);
rhs = zeros(dim*sizeBasis,1);

%assemble the stiffness matrix and RHS
elementCounter = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            elementCounter = elementCounter + 1;
            umin = PHTelem{patchIndex}(i).vertex(1);
            umax = PHTelem{patchIndex}(i).vertex(2);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (umax - umin)/2;
            
            nument = size(PHTelem{patchIndex}(i).C,1);
            sctrx = PHTelem{patchIndex}(i).nodes(1:nument); %replace with nodesGlobal for MP
            
            localstiff = zeros(nument, nument); %local stiffness
            localrhs = zeros(nument, 1); %local RHS vector
            
            %loop over the ngauss_x points on each element
            for ii=1:ngauss_x
                
                %evaluate the derivatives of the mapping from parameter
                %space to physical space
                
                [coord, dxdxi] = paramMap1D(GeoMesh{patchIndex}, gaussCoordU(ii), umin, umax);
                [a0, a1] = coefPDE(coord);

                J = det(dxdxi);
                
                %evaluate the basis functions and derivatives
                R = (PHTelem{patchIndex}(i).C)*Bu(ii,:)';
                dRdx = (PHTelem{patchIndex}(i).C)*dBdu(ii,:)';
                
                %multiply by the jacobian of the transformation from reference
                %space to the parameter space
                dRdx = dRdx*2/(umax-umin);
                
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dRdx;                                
                
                % evaluate the RHS function and assemble the local stiffness matrix                
                [~, ~, P] = exact_sol1d(coord);                                
                localstiff = localstiff +  (dR*dR'+a1*R*dR'+a0*R*R') * scalefac * gaussWeightU(ii).*J;
                localrhs = localrhs + R*P*J*gaussWeightU(ii)*scalefac;
            end
            
            %stiff(scrtx, scrtx) = stiff(scrtx, scrtx) + localstiff;
            II(indexCounter+1:indexCounter+nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter+nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter+nument^2) = reshape(localstiff,1,nument^2);
            indexCounter = indexCounter + nument^2;
            rhs(sctrx) = rhs(sctrx) + localrhs;
        end
    end
end

stiff = sparse(II,JJ,S,sizeBasis,sizeBasis);
disp(['The mesh has ', num2str(elementCounter), ' active elements.'])
