function [quadRef,estErrorGlob]=recoverDerivEstGalIsoMPDorfler(PHTelem, controlPts, sol0, quadList, p, q, Cmat, targetScale)
%recover by refining the elements with largest contribution to error, not relative contribution
%function for Galerkin method
%supports multipatches

numPatches = length(PHTelem);
estError = cell(1, numPatches); %return an array containing the estimated error in each quad
quadRef = cell(1, numPatches); %flag which quads to be refined

invC = inv(Cmat); %take the inverse of elasticity matrix (for energy norm calculations)
%c_net = PHTmesh.c_net;
numGaussU = p+1;
numGaussV = q+1;
[gaussWeightsU, gaussLocationsU] = quadrature(numGaussU, 'GAUSS', 1);
[gaussWeightsV, gaussLocationsV] = quadrature(numGaussV, 'GAUSS', 1);

err_glob = 0;
h1norm_glob = 0;
toleq = 1e-10;

up = returnUEst2alpha(p);
numQuads = 0;

for indexPatch=1:numPatches
    estError{indexPatch} = zeros(1, length(PHTelem{indexPatch}));
    quadRef{indexPatch} = zeros(1, size(quadList{indexPatch},1));
    for j=1:size(quadList{indexPatch},1)
        
        %map the Ultra-convergent points from [-1,1] to [0..1].
        uConvPts_u = ([up, up(2:end)+2]+1)/4;
        uConvPts_v = ([up, up(2:end)+2]+1)/4;
        
        knotU = [zeros(1,2*p-2),1/2,ones(1,2*p-2)];
        knotV = [zeros(1,2*q-2),1/2,ones(1,2*q-2)];
        
        p_loc = 2*p-3;
        q_loc = 2*q-3;
        
        lenU = length(knotU)-p_loc-1;  %number of basis functions in the u direction
        
        elementCounterLoc = 0;
        elementVertexLoc = zeros(4, 4);
        elementNodeLoc = zeros(4, (p_loc+1)*(q_loc+1));
        %calculate the elementNodes_loc, pt_elem
        for jj=1:length(knotV)-1
            for ii=1:length(knotU)-1
                if (abs(knotU(ii+1)-knotU(ii))>toleq) && (abs(knotV(jj+1)-knotV(jj))>toleq)
                    elementCounterLoc = elementCounterLoc + 1;
                    elementVertexLoc(elementCounterLoc, :) = [knotU(ii), knotV(jj), knotU(ii+1), knotV(jj+1)];
                    tcount = 0;
                    currow = zeros(1, (p_loc+1)*(q_loc+1));
                    %now we add the nodes from i-p...i in the u direction and
                    %j-q...j in the v direction
                    for t2=jj-q_loc:jj
                        for t1 = ii-p_loc:ii
                            tcount = tcount + 1;
                            currow(tcount) = t1+(t2-1)*lenU;
                        end
                    end
                    elementNodeLoc(elementCounterLoc,:)=currow;
                end
            end
        end
        
        num_int_pts = length(uConvPts_u)*length(uConvPts_v);
        M_int = sparse(num_int_pts, num_int_pts); %the interpolation matrix
        rhs_int = zeros(num_int_pts,3); %interpolation RHS for the stresses (xx, yy and xy components)
        
        C_loc_u = bezierExtraction(knotU, p_loc);
        C_loc_v = bezierExtraction(knotV, q_loc);
        C_loc = zeros(size(C_loc_v,1)*size(C_loc_u,1), size(C_loc_v,2)*size(C_loc_u,2), 4);
        C_loc(:,:,1) = kron(C_loc_v(:,:,1),C_loc_u(:,:,1));
        C_loc(:,:,2) = kron(C_loc_v(:,:,1),C_loc_u(:,:,2));
        C_loc(:,:,3) = kron(C_loc_v(:,:,2),C_loc_u(:,:,1));
        C_loc(:,:,4) = kron(C_loc_v(:,:,2),C_loc_u(:,:,2));
        
        row_index = 0;
        for jj=1:length(uConvPts_v)
            for ii=1:length(uConvPts_u)
                row_index = row_index + 1;
                ucoord = uConvPts_u(ii);
                vcoord = uConvPts_v(jj);
                if (vcoord >= 0) && (vcoord <= 0.5)
                    v_hat = 4*vcoord-1;
                    if(ucoord >=0) && (ucoord <= 0.5)
                        e = 1;
                        u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
                    else
                        e = 2;
                        u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
                    end
                else
                    v_hat = 4*vcoord-3;
                    if(ucoord >=0) && (ucoord <= 0.5)
                        e = 3;
                        u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
                    else
                        e = 4;
                        u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
                    end
                end
                e_glob = quadList{indexPatch}(j,e);
                
                nodes_loc = elementNodeLoc(e,:);
               
                xmin = PHTelem{indexPatch}(e_glob).vertex(1);
                xmax = PHTelem{indexPatch}(e_glob).vertex(3);
                ymin = PHTelem{indexPatch}(e_glob).vertex(2);
                ymax = PHTelem{indexPatch}(e_glob).vertex(4);
                
                N = phtBasis(u_hat, v_hat, C_loc(:,:,e), p_loc, q_loc);
                
                nument = size(PHTelem{indexPatch}(e_glob).C,1); %number of basis functions with support on current knotspan
                scrtx = PHTelem{indexPatch}(e_glob).nodesGlobal(1:nument);
                dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
                nodes = PHTelem{indexPatch}(e_glob).nodes(1:nument);
                cpts = controlPts{indexPatch}(nodes, 1:2);
                wgts = controlPts{indexPatch}(nodes, 3);
                
                Ce = PHTelem{indexPatch}(e_glob).C;
                [~, dR] = phtBasisIso(u_hat, v_hat, Ce, p, q, wgts);
                % dR
                %multiply by the jacobian of the transformation from reference
                %space to the parameter space
                dR(1,:) = dR(1,:)*2/(xmax-xmin);
                dR(2,:) = dR(2,:)*2/(ymax-ymin);
                dxdxi = dR*cpts;
                % dxdxi
                % Solve for first derivatives in global coordinates
                if abs(det(dxdxi))>1e-10
                    dR = dxdxi\dR;
                else
                    dR = pinv(dxdxi)*dR;
                end
                % dR
                B = zeros(2*nument,3);
                B(1:2:2*nument-1,1) = dR(1,:)';
                B(2:2:2*nument,2) = dR(2,:)';
                B(1:2:2*nument-1,3) = dR(2,:)';
                B(2:2:2*nument,3) = dR(1,:)';
                
                stressvect = Cmat*B'*sol0(dscrtx);
                
                %set-up linear system for solving for the recovered stress
                M_int(row_index,nodes_loc) = N;
                rhs_int(row_index,:) = stressvect';
                
            end
        end
        %M_int
        %rhs_int
        rec_sol = M_int\rhs_int;
        %pause
        %calculate the error in energy norm for each element of the patch
        for e=1:4 %for each element in the patch
            
            e_glob = quadList{indexPatch}(j,e);
            
            nument = size(PHTelem{indexPatch}(e_glob).C,1); %number of basis functions with support on current knotspan
            scrtx = PHTelem{indexPatch}(e_glob).nodesGlobal(1:nument);
            dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2*nument);
            nodes = PHTelem{indexPatch}(e_glob).nodes(1:nument);
            cpts = controlPts{indexPatch}( nodes, 1:2);
            wgts = controlPts{indexPatch}(nodes, 3);
            nodes_loc = elementNodeLoc(e,:);
            
            Ce = PHTelem{indexPatch}(e_glob).C;
            Ce_loc = C_loc(:,:,e);
            
            xmin = PHTelem{indexPatch}(e_glob).vertex(1);
            xmax = PHTelem{indexPatch}(e_glob).vertex(3);
            ymin = PHTelem{indexPatch}(e_glob).vertex(2);
            ymax = PHTelem{indexPatch}(e_glob).vertex(4);
            
            detJacobian = (xmax - xmin)*(ymax - ymin)/4;
            
            h1normerr = 0;
            for jj=1:numGaussV
                for ii=1:numGaussU
                    u_hat = gaussLocationsU(ii); %gauss points on [-1, 1]
                    v_hat = gaussLocationsV(jj);
                    
                    [~, dN] = phtBasisIso(u_hat, v_hat, Ce, p, q, wgts);
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dN(1,:) = dN(1,:)*2/(xmax-xmin);
                    dN(2,:) = dN(2,:)*2/(ymax-ymin);
                    
                    dxdxi = dN*cpts;
                    J = det(dxdxi);
                    
                    % Solve for first derivatives in global coordinates
                    dN = dxdxi\dN;
                    
                    shape = phtBasis(u_hat, v_hat, Ce_loc, p_loc, q_loc);
                    recovStress = shape*rec_sol(nodes_loc,:);
                    
                    B = zeros(2*nument,3);
                    B(1:2:2*nument-1,1) = dN(1,:)';
                    B(2:2:2*nument,2) = dN(2,:)';
                    B(1:2:2*nument-1,3) = dN(2,:)';
                    B(2:2:2*nument,3) = dN(1,:)';
                    % sol0(dscrtx)
                    compStress = Cmat*B'*sol0(dscrtx);
                    recovStress = recovStress';
                    h1norm_glob = h1norm_glob + recovStress'*invC*recovStress*gaussWeightsU(ii)*gaussWeightsV(jj)*detJacobian*J;
                    h1normerr = h1normerr + (recovStress-compStress)'*invC*(recovStress-compStress)*gaussWeightsU(ii)*gaussWeightsV(jj)*detJacobian*J;
                end
            end
            err_glob = err_glob + h1normerr;
            estError{indexPatch}(e_glob) = h1normerr;
        end
    end
    numQuads = numQuads+size(quadList{indexPatch},1);
end
%estError

estErrorGlob = sqrt(err_glob/h1norm_glob);
quadCounter = 0;
quadErrors = zeros(1, numQuads);
quadIndex = zeros(numQuads,2);
for indexPatch = 1:numPatches
    for j=1:size(quadList{indexPatch},1)
        quadCounter = quadCounter + 1;
        quadErrors(quadCounter) = sum(estError{indexPatch}(quadList{indexPatch}(j,:)));        
        quadIndex(quadCounter,:) = [indexPatch, j];
    end
end

[quadErrors, index] = sort(quadErrors,'descend');
quadIndex = quadIndex(index, :);
%pause
targetError = targetScale*sqrt(err_glob);
sumError = 0;
for i=1:numQuads
    sumError = sumError + quadErrors(i);
    quadRef{quadIndex(i,1)}(quadIndex(i,2)) = 1;
    if sqrt(sumError) > targetError
        disp(['Refining ', num2str(i), ' quads with largest error out of ', num2str(numQuads)])
        break
    end
end
