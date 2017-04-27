function [ stiff, mass, rhs ] = assembleGalerkinMassSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, r, Cmat, modelRho, octupleList )
%assembles the stiffness matrix and rhs (Galerkin method)
%uses GIFT mapping
%supports multipatches


%Gauss points
ngauss_x = p+1;
ngauss_y = q+1;
ngauss_z = r+1;
[gauss_weight_x, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[gauss_weight_y, gauss_coord_y] = quadrature( ngauss_y, 'GAUSS', 1 );
[gauss_weight_z, gauss_coord_z] = quadrature( ngauss_y, 'GAUSS', 1 );


%take the transpose so that they are in the format expected by
%bernstein_basis
gauss_coord_x = gauss_coord_x';
gauss_coord_y = gauss_coord_y';
gauss_coord_z = gauss_coord_z';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v] = bernstein_basis(gauss_coord_y,q);
[B_w, dB_w] = bernstein_basis(gauss_coord_z,r);

R = zeros(ngauss_x, ngauss_y, ngauss_z, (p+1)*(q+1)*(r+1));
dBdu = zeros(ngauss_x, ngauss_y, ngauss_z, (p+1)*(q+1)*(r+1));
dBdv = zeros(ngauss_x, ngauss_y, ngauss_z, (p+1)*(q+1)*(r+1));
dBdw = zeros(ngauss_x, ngauss_y, ngauss_z, (p+1)*(q+1)*(r+1));

%the derivatives of the 3D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:ngauss_z
                for jj=1:ngauss_y
                    for ii=1:ngauss_x
                        R(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdu(ii,jj,kk,basisCounter) = dB_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdv(ii,jj,kk,basisCounter) = B_u(ii,i)*dB_v(jj,j)*B_w(kk,k);
                        dBdw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*dB_w(kk,k);
                    end
                end
            end
        end
    end
end


%initialize LHS stiffness matrix and RHS vector

dim = 3; %the dimension of physical space
stiff = sparse(dim*sizeBasis,dim*sizeBasis);
mass = sparse(dim*sizeBasis,dim*sizeBasis);
rhs = zeros(dim*sizeBasis,1);
%assemble the collocation matrix and RHS

%ass_time_counter = 0;
for patchIndex = 1:length(PHTelem)
    disp(['Assembling patch ', num2str(patchIndex), ' of ', num2str(length(PHTelem)), '...'])
    for octIndex = 1:size(octupleList{patchIndex},1)
        % toc
        disp(['Assembling octuple ', num2str(octIndex), ' of ', num2str(size(octupleList{patchIndex},1)),'...'])
        nument = (p+1)*(q+1)*(r+1);
        element_block_size = 8;
        II = zeros(1, element_block_size*nument^2*dim^2);
        JJ = zeros(1, element_block_size*nument^2*dim^2);
        indexcounter = 0;
        S = zeros(1, element_block_size*nument^2*dim^2);
        
        for i=octupleList{patchIndex}(octIndex,:)
            
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(4);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(5);
            zmin = PHTelem{patchIndex}(i).vertex(3);
            zmax = PHTelem{patchIndex}(i).vertex(6);
            
            %the jacobian of the transformation from [-1,1]x[-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax] x [zmin, zmax]
            scalefac = (xmax - xmin)*(ymax - ymin)*(zmax-zmin)/8;
            
            nument = size(PHTelem{patchIndex}(i).C,1);
            scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
            B = zeros(6, dim*nument);
            
            localstiff = zeros(dim*nument, dim*nument); %local stiffness
            localmass = zeros(dim*nument, dim*nument); %local mass
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            
            
            for kk=1:ngauss_z
                for jj=1:ngauss_y
                    for ii=1:ngauss_x
                        
                        %evaluate the derivatives of the mapping from parameter
                        %space to physical space
                        [~, dxdxi] = paramMap3D( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), gauss_coord_z(kk), xmin, ymin, zmin, xmax, ymax, zmax);
                        %
                        J = abs(det(dxdxi));
                        
                        cR = (PHTelem{patchIndex}(i).C)*squeeze(R(ii,jj,kk,:));
                        dRdx = (PHTelem{patchIndex}(i).C)*squeeze(dBdu(ii,jj,kk,:));
                        dRdy = (PHTelem{patchIndex}(i).C)*squeeze(dBdv(ii,jj,kk,:));
                        dRdz = (PHTelem{patchIndex}(i).C)*squeeze(dBdw(ii,jj,kk,:));
                        
                        %multiply by the jacobian of the transformation from reference
                        %space to the parameter space
                        dRdx = dRdx*2/(xmax-xmin);
                        dRdy = dRdy*2/(ymax-ymin);
                        dRdz = dRdz*2/(zmax-zmin);
                        
                        % Solve for first derivatives in global coordinates
                        dR = dxdxi\[dRdx';dRdy';dRdz'];
                        
                        B(1,1:3:3*nument-2) = dR(1,:);
                        B(2,2:3:3*nument-1) = dR(2,:);
                        B(3,3:3:3*nument) = dR(3,:);
                        
                        B(4,1:3:3*nument-2) = dR(2,:);
                        B(4,2:3:3*nument-1) = dR(1,:);
                        
                        B(5,2:3:3*nument-1) = dR(3,:);
                        B(5,3:3:3*nument) = dR(2,:);
                        
                        B(6,1:3:3*nument-2) = dR(3,:);
                        B(6,3:3:3*nument) = dR(1,:);
                        
                        % N matrix
                        N(1,1:3:3*nument-2) = cR;
                        N(2,2:3:3*nument-1) = cR;
                        N(3,3:3:3*nument)   = cR;
                        
                        %TODO: implement non-zero volume force
                        localstiff = localstiff + B' * Cmat * B * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj)*gauss_weight_z(kk).*J;
                        localmass = localmass +  N' * (modelRho * N) * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj)*gauss_weight_z(kk).*J;
                    end
                end
            end
            
            II(indexcounter+1:indexcounter+dim^2*nument^2) = repmat(tscrtx,1,dim*nument);
            JJ(indexcounter+1:indexcounter+dim^2*nument^2) = reshape(repmat(tscrtx',1,dim*nument)',1,dim^2*nument^2);
            S(indexcounter+1:indexcounter+dim^2*nument^2) = reshape(localstiff,1,dim^2*nument^2);
            SM(indexcounter+1:indexcounter+dim^2*nument^2) = reshape(localmass,1,dim^2*nument^2);
            indexcounter = indexcounter + dim^2*nument^2;
            
        end
        stiff = stiff + sparse(II,JJ,S,dim*sizeBasis,dim*sizeBasis);
        mass = mass + sparse(II,JJ,SM,dim*sizeBasis,dim*sizeBasis);
    end
end