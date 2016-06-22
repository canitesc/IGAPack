function plotStressDisp3DMP(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
%plots the computed solution
%supports multipatches

numPlotX = 2;
numPlotY = 2;
numPlotZ = 2;

edgeSpace = 1e-3; %space between the plot cell and the true element boundary
noGpEle = numPlotX*numPlotY*numPlotZ;

%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            numElem = numElem+1;
        end
    end
end


noElems =  numElem;
x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % displacements of Gauss points
sigma = zeros(6,noGpEle,noElems);  % stresses      of Gauss points

plotPtsX = [-1+edgeSpace,1-edgeSpace];
plotPtsY = [-1+edgeSpace,1-edgeSpace];
plotPtsZ = [-1+edgeSpace,1-edgeSpace];

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(plotPtsX,p);
[B_v, dB_v] = bernstein_basis(plotPtsY,q);
[B_w, dB_w] = bernstein_basis(plotPtsZ,r);

B_uvw = zeros(numPlotX, numPlotY, numPlotZ, (p+1)*(q+1)*(r+1));
dBdu = zeros(numPlotX, numPlotY, numPlotZ, (p+1)*(q+1)*(r+1));
dBdv = zeros(numPlotX, numPlotY, numPlotZ, (p+1)*(q+1)*(r+1));
dBdw = zeros(numPlotX, numPlotY, numPlotZ, (p+1)*(q+1)*(r+1));

%compute the function values and derivatives of the 3D Bernstein polynomials at plot points on the
%master element
basisCounter = 0;
for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            for kk=1:numPlotZ
                for jj=1:numPlotY
                    for ii=1:numPlotX
                        B_uvw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdu(ii,jj,kk,basisCounter) = dB_u(ii,i)*B_v(jj,j)*B_w(kk,k);
                        dBdv(ii,jj,kk,basisCounter) = B_u(ii,i)*dB_v(jj,j)*B_w(kk,k);
                        dBdw(ii,jj,kk,basisCounter) = B_u(ii,i)*B_v(jj,j)*dB_w(kk,k);
                    end
                end
            end
        end
    end
end

e = 1;
for patchIndex = 1:length(PHTelem)
    
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(4);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(5);
            zmin = PHTelem{patchIndex}(i).vertex(3);
            zmax = PHTelem{patchIndex}(i).vertex(6);
            
            nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan
            scrt = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            
            scrt_x = 3*scrt-2;
            scrt_y = 3*scrt-1;
            scrt_z = 3*scrt;
            tscrtx = reshape([scrt_x; scrt_y; scrt_z],1,3*nument);
            
            
            gp = 1;
            
            B = zeros(6, 3*nument);
            dR  = zeros(3, nument);
            
            for kk=1:numPlotZ
                for jj=1:numPlotY
                    for ii=1:numPlotX
                        
                        [coord, dxdxi] = paramMap3D( GIFTmesh{patchIndex}, plotPtsX(ii), plotPtsY(jj), plotPtsZ(kk), xmin, ymin, zmin, xmax, ymax, zmax);
                        
                        cR = PHTelem{patchIndex}(i).C * squeeze(B_uvw(ii,jj,kk,:));
                        
                        dR(1,:) = PHTelem{patchIndex}(i).C * squeeze(dBdu(ii,jj,kk,:))*2/(xmax-xmin);
                        dR(2,:) = PHTelem{patchIndex}(i).C * squeeze(dBdv(ii,jj,kk,:))*2/(ymax-ymin);
                        dR(3,:) = PHTelem{patchIndex}(i).C * squeeze(dBdw(ii,jj,kk,:))*2/(zmax-zmin);
                        
                        % Solve for first derivatives in global coordinates
                        dR = pinv(dxdxi)*dR;
                        
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
                        disp_x = cR'*sol0(scrt_x);
                        disp_y = cR'*sol0(scrt_y);
                        disp_z = cR'*sol0(scrt_z);
                        
                        stressvect = Cmat*B*sol0(tscrtx);
                        
                        x(1,gp,e)     = coord(1);
                        x(2,gp,e)     = coord(2);
                        x(3,gp,e)     = coord(3);
                        
                        u(1,gp,e)     = disp_x;
                        u(2,gp,e)     = disp_y;
                        u(3,gp,e)     = disp_z;
                        sigma(1:6,gp,e) = stressvect;
                        gp = gp +1 ;
                    end
                end
            end
            e = e + 1;
            
        end
    end
end
msh_to_vtu_3d(x, u, sigma, [numPlotX numPlotY numPlotZ], vtuFile);