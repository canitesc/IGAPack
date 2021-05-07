function  plotSolPHTElasticMP( sol0, PHTelem, GIFTmesh, p, q, Cmat, factor )
%supports multipatches

numPts = 2;

%plots the deformed shape + stresses
xi = linspace(-1,1,numPts);
eta = linspace(-1,1,numPts);


if nargin<7
    %set a default magnification factor for the displacement
    factor=10;
end


%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
numPatches = length(PHTelem);
for patchIndex=1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            numElem = numElem+1;
        end
    end
end


%Displaying the displacements
element4 = zeros(numElem, 4);
physcoord = zeros(4*numElem, 2);
dispcoord = zeros(4*numElem, 2);
straincoord = zeros(4*numElem, 3);
sigmacoord = zeros(4*numElem, 3);


%define the 2D Bernstein polynomials
[B_u, dB_u] = bernstein_basis(xi,p);
[B_v, dB_v] = bernstein_basis(eta,q);

Buv = zeros(numPts, numPts, numPts*numPts);
dBdu = zeros(numPts, numPts, numPts*numPts);
dBdv = zeros(numPts, numPts, numPts*numPts);

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

elementCounter = 0;
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            elementCounter = elementCounter + 1;
            element4(elementCounter, :) = (elementCounter-1)*4+1:(elementCounter-1)*4+4;
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            coord = cell(numPts,numPts);
            
            nument = size(PHTelem{indexPatch}(i).C,1); %number of basis functions with support on current knotspan
            
            %initialize matrices to store the displacement, strain and stress
            %values at each plot point
            dispmatx = zeros(numPts,numPts);
            dispmaty = zeros(numPts,numPts);
            strain11 = zeros(numPts,numPts);
            strain12 = zeros(numPts,numPts);
            strain22 = zeros(numPts,numPts);
            stressvect = cell(numPts,numPts);
            
            scrt = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            scrt_x = 2*scrt-1;
            scrt_y = 2*scrt;
            dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);
            
            for jj=1:numPts
                for ii=1:numPts
                    
                    %compute the mapping from parameter space to physical space
                    [ coord_pt, dxdxi]  = paramMap( GIFTmesh{indexPatch}, xi(ii), eta(jj), xmin, ymin, xmax, ymax);
                    
                    %evaluate the basis functions
                    
                    cR = PHTelem{indexPatch}(i).C * squeeze(Buv(ii,jj,:));
                    cdRdx = PHTelem{indexPatch}(i).C * squeeze(dBdu(ii,jj,:))*2/(xmax-xmin);
                    cdRdy = PHTelem{indexPatch}(i).C * squeeze(dBdv(ii,jj,:))*2/(ymax-ymin);
                    
                    coord{jj,ii} = [coord_pt(1), coord_pt(2)];
                    
                    % Solve for first derivatives in global coordinates
                    
                    dR = pinv(dxdxi)*[cdRdx';cdRdy'];
                    %dR = dxdxi\[cdRdx';cdRdy'];
                    
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dR(1,:);
                    B(2:2:2*nument,2) = dR(2,:);
                    B(1:2:2*nument-1,3) = dR(2,:);
                    B(2:2:2*nument,3) = dR(1,:);
                    
                    %calculate displacement values
                    dispmatx(jj,ii) = dispmatx(jj,ii) + cR'*sol0(scrt_x);
                    dispmaty(jj,ii) = dispmaty(jj,ii) + cR'*sol0(scrt_y);
                    
                    %calculate the stress values
                    stressvect{jj,ii} = Cmat*B'*sol0(dscrtx);
                    
                    
                end
            end
            
            physcoord((elementCounter-1)*4+1, :) = coord{1,1};
            physcoord((elementCounter-1)*4+2, :) = coord{1,end};
            physcoord((elementCounter-1)*4+3, :) = coord{end,end};
            physcoord((elementCounter-1)*4+4, :) = coord{end,1};
            
            dispcoord((elementCounter-1)*4+1, :) = [dispmatx(1,1) dispmaty(1,1)];
            dispcoord((elementCounter-1)*4+2, :) = [dispmatx(1,2) dispmaty(1,2)];
            dispcoord((elementCounter-1)*4+3, :) = [dispmatx(2,2) dispmaty(2,2)];
            dispcoord((elementCounter-1)*4+4, :) = [dispmatx(2,1) dispmaty(2,1)];
            straincoord((elementCounter-1)*4+1, :) = [strain11(1,1) strain12(1,1) strain22(1,1)];
            straincoord((elementCounter-1)*4+2, :) = [strain11(1,2) strain12(1,2) strain22(1,2)];
            straincoord((elementCounter-1)*4+3, :) = [strain11(2,2) strain12(2,2) strain22(2,2)];
            straincoord((elementCounter-1)*4+4, :) = [strain11(2,1) strain12(2,1) strain22(2,1)];
            
            
            sigmacoord((elementCounter-1)*4+1, :) = stressvect{1,1}';
            sigmacoord((elementCounter-1)*4+2, :) = stressvect{1,2}';
            sigmacoord((elementCounter-1)*4+3, :) = stressvect{2,2}';
            sigmacoord((elementCounter-1)*4+4, :) = stressvect{2,1}';
            
            hold on
        end
    end
end

factor = 1;

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,1), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{11}')
colorbar('vert')
drawnow

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,3), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{12}')
colorbar('vert')
drawnow

figure
trisurf(element4,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), sigmacoord(:,2), 'EdgeColor','none','facecolor','interp')
view(0,90)
title('Displacements and \sigma_{22}')
colorbar('vert')
drawnow


