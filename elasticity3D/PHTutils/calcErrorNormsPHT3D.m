function  [l2relerr, h1relerr]=calcErrorNormsPHT3D( sol0, PHTelem, GIFTmesh, p, q,r, Cmat )
%calculate the actual error in the computed solution for plate with hole
%problem
%supports multipatches

numPatches = length(PHTelem);

numGaussX = p+1;
numGaussY = q+1;
numGaussZ = r+1;

[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
[gwy, gpy]=quadrature(numGaussY, 'GAUSS', 1);
[gwz, gpz]=quadrature(numGaussZ, 'GAUSS', 1);

gpx=gpx';
gpy=gpy';
gpz=gpz';

l2norm = 0;
h1norm = 0;

l2relerr = 0;
h1relerr = 0;

invC = inv(Cmat);


%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(gpx,p);
[B_v, dB_v] = bernstein_basis(gpy,q);
[B_w, dB_w] = bernstein_basis(gpz,r);

dim = 3;

B_uv = zeros(numGaussX, numGaussY,numGaussZ,(p+1)*(q+1)*(r+1));
dBdu = zeros(numGaussX, numGaussY, numGaussZ, (p+1)*(q+1)*(r+1));
dBdv = zeros(numGaussX, numGaussY, numGaussZ, (p+1)*(q+1)*(r+1));
dBdw = zeros(numGaussX, numGaussY, numGaussZ, (p+1)*(q+1)*(r+1));

%the derivatives of the 3D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:numGaussZ
                for jj=1:numGaussY
                    for ii=1:numGaussX
                         B_uv(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdu(ii,jj,kk,basisCounter) = dB_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdv(ii,jj,kk,basisCounter) = B_u(ii,i)*dB_v(jj,j)*B_w(kk,k);
                        dBdw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*dB_w(kk,k);
                    end
                end
            end
        end
    end
end


for i=1:length(PHTelem)
    if isempty(PHTelem(i).children)
        xmin = PHTelem(i).vertex(1);
        xmax = PHTelem(i).vertex(4);
        ymin = PHTelem(i).vertex(2);
        ymax = PHTelem(i).vertex(5);
        zmin = PHTelem(i).vertex(3);
        zmax = PHTelem(i).vertex(6);
            
           %the jacobian of the transformation from [-1,1]x[-1,1]x[-1,1] to
        %[xmin, xmax]x [ymin, ymax] x [zmin, zmax]
        scalefac = (xmax - xmin)*(ymax - ymin)*(zmax-zmin)/8;
            
            
           % J = zeros(numGaussY, numGaussX, numGaussZ);
            
            
            nument = size(PHTelem(i).C,1);
            scrt = PHTelem(i).nodes(1:nument);
            scrtx = 3*scrt-2;
            scrty = 3*scrt-1;
            scrtz = 3*scrt;
            tscrtx = reshape([3*scrt-2; 3*scrt-1; 3*scrt],1,3*nument);
            B = zeros(6, dim*nument);
            
             %loop over the ngauss_x x ngauss_y gauss points on each element
        
        
        for kk=1:numGaussZ
            for jj=1:numGaussY
                for ii=1:numGaussX
                    
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    [coord, dxdxi] = paramMap3D( GIFTmesh, gpx(ii), gpy(jj), gpz(kk), xmin, ymin, zmin, xmax, ymax, zmax);                     
                    %
                   [ disp_ex, stress_ex ] = exactSolHemi( coord );
                    
                    J = det(dxdxi);
                    dRdx = (PHTelem(i).C)*squeeze(dBdu(ii,jj,kk,:));
                    dRdy = (PHTelem(i).C)*squeeze(dBdv(ii,jj,kk,:));
                    dRdz = (PHTelem(i).C)*squeeze(dBdw(ii,jj,kk,:));
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    dRdz = dRdz*2/(zmax-zmin);
                    
                     cR = (PHTelem(i).C)* squeeze(B_uv(ii,jj,kk,:));
                    
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
                    
                     %calculate displacement values
                    disp_x = cR'*sol0(scrtx);
                    disp_y = cR'*sol0(scrty);
                    disp_z = cR'*sol0(scrtz);
                    
                    %calculate the error in stress values
                    stressvect = Cmat*B*sol0(tscrtx);
               %     stress_ex
               %     stressvect'
               %     pause
                    
                    l2norm = l2norm + (disp_ex(1)^2 + disp_ex(2)^2 + disp_ex(3)^2)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                    h1norm = h1norm + stress_ex*invC*stress_ex'*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                    
                    l2relerr = l2relerr + ((disp_ex(1)-disp_x)^2 + (disp_ex(2)-disp_y)^2 + (disp_ex(3)-disp_z)^2)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                    h1relerr = h1relerr + (stress_ex-stressvect')*invC*(stress_ex'-stressvect)*gwx(ii)*gwy(jj)*gwz(kk)*scalefac*J;
                end
            end
        end
    end
end

l2relerr = sqrt(l2relerr)/sqrt(l2norm);
h1relerr = sqrt(h1relerr)/sqrt(h1norm);


            
%             for jj=1:numGaussY
%                 for ii=1:numGaussX
%                     
%                     [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
%                     
%                     disp_ex = holeu_d([coord(1), coord(2)], rad, Emod, nu, tx);
%                     stress_ex = ghole(coord(1), coord(2), rad, tx);
%                     
%                     J(jj, ii) = det(dxdxi);
%                     cR = PHTelem{patchIndex}(i).C * squeeze(B_uv(ii,jj,:));
%                     
%                     dR(1,:) = PHTelem{patchIndex}(i).C * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
%                     dR(2,:) = PHTelem{patchIndex}(i).C * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
%                     
%                     % Solve for first derivatives in global coordinates
%                     dR = dxdxi\dR;
%                     
%                     B = zeros(2*nument,3);
%                     B(1:2:2*nument-1,1) = dR(1,:);
%                     B(2:2:2*nument,2) = dR(2,:);
%                     B(1:2:2*nument-1,3) = dR(2,:);
%                     B(2:2:2*nument,3) = dR(1,:);
%                     
%                     %calculate displacement values
%                     disp_x = cR'*sol0(scrt_x);
%                     disp_y = cR'*sol0(scrt_y);
%                     
%                     %calculate the error in stress values
%                     stressvect = Cmat*B'*sol0(dscrtx);
%                     
%                     l2norm = l2norm + (disp_ex(1)^2 + disp_ex(2)^2)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
%                     h1norm = h1norm + stress_ex'*invC*stress_ex*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
%                     
%                     l2relerr = l2relerr + ((disp_ex(1)-disp_x)^2 + (disp_ex(2)-disp_y)^2)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
%                     h1relerr = h1relerr + (stress_ex'-stressvect')*invC*(stress_ex-stressvect)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
%                     
%                 end
%             end
%             
%         end
% end
% 
% l2relerr = sqrt(l2relerr)/sqrt(l2norm);
% h1relerr = sqrt(h1relerr)/sqrt(h1norm);
