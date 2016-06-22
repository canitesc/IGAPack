function  [l2relerr, h1relerr]=calcErrorNormsPH( sol0, PHTelem, GIFTmesh, p, q, Cmat, Emod, nu, rad, tx )
%calculate the actual error in the computed solution for plate with hole
%problem
%supports multipatches

numPatches = length(PHTelem);

numGaussX = p+1;
numGaussY = q+1;

[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
[gwy, gpy]=quadrature(numGaussY, 'GAUSS', 1);

gpx=gpx';
gpy=gpy';

l2norm = 0;
h1norm = 0;

l2relerr = 0;
h1relerr = 0;

invC = inv(Cmat);


%define the 2D Bernstein polynomials
[B_u, dB_u] = bernstein_basis(gpx,p);
[B_v, dB_v] = bernstein_basis(gpy,q);
B_uv = zeros(numGaussX, numGaussY,(p+1)*(q+1));
dBdu = zeros(numGaussX, numGaussY,(p+1)*(q+1));
dBdv = zeros(numGaussX, numGaussY,(p+1)*(q+1));

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        B_uv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            %jacobian of the transformation from reference [-1,1]x[-1,1]
            %element to the local element in parameter
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            
            
            J = zeros(numGaussY, numGaussX);
            
            
            nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan
            
            scrt = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            scrt_x = 2*scrt-1;
            scrt_y = 2*scrt;
            dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);
            
            for jj=1:numGaussY
                for ii=1:numGaussX
                    
                    [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
                    
                    disp_ex = holeu_d([coord(1), coord(2)], rad, Emod, nu, tx);
                    stress_ex = ghole(coord(1), coord(2), rad, tx);
                    
                    J(jj, ii) = det(dxdxi);
                    cR = PHTelem{patchIndex}(i).C * squeeze(B_uv(ii,jj,:));
                    
                    dR(1,:) = PHTelem{patchIndex}(i).C * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                    dR(2,:) = PHTelem{patchIndex}(i).C * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                    
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\dR;
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dR(1,:);
                    B(2:2:2*nument,2) = dR(2,:);
                    B(1:2:2*nument-1,3) = dR(2,:);
                    B(2:2:2*nument,3) = dR(1,:);
                    
                    %calculate displacement values
                    disp_x = cR'*sol0(scrt_x);
                    disp_y = cR'*sol0(scrt_y);
                    
                    %calculate the error in stress values
                    stressvect = Cmat*B'*sol0(dscrtx);
                    
                    l2norm = l2norm + (disp_ex(1)^2 + disp_ex(2)^2)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
                    h1norm = h1norm + stress_ex'*invC*stress_ex*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
                    
                    l2relerr = l2relerr + ((disp_ex(1)-disp_x)^2 + (disp_ex(2)-disp_y)^2)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
                    h1relerr = h1relerr + (stress_ex'-stressvect')*invC*(stress_ex-stressvect)*gwx(ii)*gwy(jj)*scalefac*J(jj,ii);
                    
                end
            end
            
        end
    end
end
l2relerr = sqrt(l2relerr)/sqrt(l2norm);
h1relerr = sqrt(h1relerr)/sqrt(h1norm);
