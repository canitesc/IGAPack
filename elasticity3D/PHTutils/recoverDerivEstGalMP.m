function [octupleRef,estErrorGlob, estError]=recoverDerivEstGalMP(PHTelem, GIFTmesh, sol0, target_rel_error, octupleList, p, q, r, Cmat, targetScale)
%recover by refining the elements with largest contribution to error, not relative contribution
%function for Galerkin method
%assumes p=q=r

alpha = floor((p-1)/2); 


estError = zeros(1, length(PHTelem)); %return an array containing the estimated error in each patch
octupleRef = zeros(1, size(octupleList,1)); %flag which patches to be refined

invC = inv(Cmat); %take the inverse of elasticity matrix (for energy norm calculations)
%c_net = PHTmesh.c_net;
numGaussU = p+1;
numGaussV = q+1;
numGaussW = r+1;
[gaussWeightsU, gaussLocationsU] = quadrature(numGaussU, 'GAUSS', 1);
[gaussWeightsV, gaussLocationsV] = quadrature(numGaussV, 'GAUSS', 1);
[gaussWeightsW, gaussLocationsW] = quadrature(numGaussW, 'GAUSS', 1);

err_glob = 0;
h1norm_glob = 0;
toleq = 1e-10;

up = returnUEst(p);

%map the Ultra-convergent points from [-1,1] to [0..1].

if floor(p/2)==(p/2)
    uConvPts_u = ([up, up+2]+1)/4;
    uConvPts_v = ([up, up+2]+1)/4;
    uConvPts_w = ([up, up+2]+1)/4;
else
    uConvPts_u = ([up, up(2:end)+2]+1)/4;
    uConvPts_v = ([up, up(2:end)+2]+1)/4;
    uConvPts_w = ([up, up(2:end)+2]+1)/4;
end
p_loc = p;
q_loc = q;
r_loc = r;

knotU = [zeros(1,p_loc+1),1/2*ones(1,p-alpha-1),ones(1,p_loc+1)];
knotV = [zeros(1,q_loc+1),1/2*ones(1,q-alpha-1),ones(1,q_loc+1)];
knotW = [zeros(1,r_loc+1),1/2*ones(1,r-alpha-1),ones(1,r_loc+1)];

lenU = length(knotU)-p_loc-1;  %number of basis functions in the u direction
lenV = length(knotV)-q_loc-1;  %number of basis functions in the v direction
lenW = length(knotW)-r_loc-1;  %number of basis functions in the w direction

elementCounterLoc = 0;
elementNodeLoc = zeros(4, (p_loc+1)*(q_loc+1)*(r_loc+1));
%calculate the elementNodes_loc, pt_elem
for kk=1:length(knotW)-1
    for jj=1:length(knotV)-1
        for ii=1:length(knotU)-1
            if (abs(knotU(ii+1)-knotU(ii))>toleq) && (abs(knotV(jj+1)-knotV(jj))>toleq) && (abs(knotW(kk+1)-knotW(kk))>toleq)
                elementCounterLoc = elementCounterLoc + 1;
                
                tcount = 0;
                currow = zeros(1, (p_loc+1)*(q_loc+1)*(r_loc+1));
                %now we add the nodes from i-p...i in the u
                %direction, j-q...j in the v direction, k-r...k in
                %w direction
                for t3 = kk-r_loc:kk
                    for t2=jj-q_loc:jj
                        for t1 = ii-p_loc:ii
                            tcount = tcount + 1;
                            currow(tcount) = t1+(t2-1)*lenU+(t3-1)*lenU*lenV;
                        end
                    end
                end
                elementNodeLoc(elementCounterLoc,:)=currow;
            end
        end
    end
end

C_loc_u = bezierExtraction(knotU, p_loc);
C_loc_v = bezierExtraction(knotV, q_loc);
C_loc_w = bezierExtraction(knotW, r_loc);
%C_loc = zeros(size(C_loc_w,1)*size(C_loc_v,1)*size(C_loc_u,1), size(C_loc_w,2)*size(C_loc_v,2)*size(C_loc_u,2), 8);


C_loc(:,:,1) =  kron(kron(C_loc_w(:,:,1),C_loc_v(:,:,1)),C_loc_u(:,:,1));
C_loc(:,:,2) =  kron(kron(C_loc_w(:,:,1),C_loc_v(:,:,1)),C_loc_u(:,:,2));
C_loc(:,:,3) =  kron(kron(C_loc_w(:,:,1),C_loc_v(:,:,2)),C_loc_u(:,:,1));
C_loc(:,:,4) =  kron(kron(C_loc_w(:,:,1),C_loc_v(:,:,2)),C_loc_u(:,:,2));
C_loc(:,:,5) =  kron(kron(C_loc_w(:,:,2),C_loc_v(:,:,1)),C_loc_u(:,:,1));
C_loc(:,:,6) =  kron(kron(C_loc_w(:,:,2),C_loc_v(:,:,1)),C_loc_u(:,:,2));
C_loc(:,:,7) =  kron(kron(C_loc_w(:,:,2),C_loc_v(:,:,2)),C_loc_u(:,:,1));
C_loc(:,:,8) =  kron(kron(C_loc_w(:,:,2),C_loc_v(:,:,2)),C_loc_u(:,:,2));



for j=1:size(octupleList,1)
    disp(['Estimating the error in octuple ', num2str(j), ' of ', num2str(size(octupleList,1)), '...'])
    
    num_int_pts = length(uConvPts_u)*length(uConvPts_v)*length(uConvPts_w);
    M_int = sparse(num_int_pts, lenU*lenV*lenW); %the interpolation matrix
    rhs_int = zeros(num_int_pts,6); %interpolation RHS for the stresses (xx, yy, zz, xy, yz, xz components)
    
    row_index = 0;
    for kk=1:length(uConvPts_w)
        for jj=1:length(uConvPts_v)
            for ii=1:length(uConvPts_u)
                row_index = row_index + 1;
                ucoord = uConvPts_u(ii);
                vcoord = uConvPts_v(jj);
                wcoord = uConvPts_w(kk);
                if (wcoord >= 0) && (wcoord <= 0.5)
                    w_hat = 4*wcoord-1;
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
                else
                    w_hat = 4*wcoord-3;
                    if (vcoord >= 0) && (vcoord <= 0.5)
                        v_hat = 4*vcoord-1;
                        if(ucoord >=0) && (ucoord <= 0.5)
                            e = 5;
                            u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
                        else
                            e = 6;
                            u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
                        end
                    else
                        v_hat = 4*vcoord-3;
                        if(ucoord >=0) && (ucoord <= 0.5)
                            e = 7;
                            u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
                        else
                            e = 8;
                            u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
                        end
                    end
                    
                end
                e_glob = octupleList(j,e);
                
                nodes_loc = elementNodeLoc(e,:);
                
                xmin = PHTelem(e_glob).vertex(1);
                xmax = PHTelem(e_glob).vertex(4);
                ymin = PHTelem(e_glob).vertex(2);
                ymax = PHTelem(e_glob).vertex(5);
                zmin = PHTelem(e_glob).vertex(3);
                zmax = PHTelem(e_glob).vertex(6);
                
                N = phtBasis(u_hat, v_hat, w_hat, C_loc(:,:,e), p_loc, q_loc, r_loc);
                
                Ce = PHTelem(e_glob).C;
                
                [ ~, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
                [~, dR] = phtBasis(u_hat, v_hat, w_hat, Ce, p, q, r);
                
                
                nument = size(PHTelem(e_glob).C,1); %number of basis functions with support on current knotspan
                scrtx = PHTelem(e_glob).nodesGlobal(1:nument);
                tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
                
                %multiply by the jacobian of the transformation from reference
                %space to the parameter space
                dR(1,:) = dR(1,:)*2/(xmax-xmin);
                dR(2,:) = dR(2,:)*2/(ymax-ymin);
                dR(3,:) = dR(3,:)*2/(zmax-zmin);
                
                % Solve for first derivatives in global coordinates
                if abs(det(dxdxi))>1e-10
                    dR = dxdxi\dR;
                else
                    dR = pinv(dxdxi)*dR;
                end
                
                %initialize the B-matrix for 3D
                B = zeros(6,3*nument);
                
                B(1,1:3:3*nument-2) = dR(1,:)';
                B(2,2:3:3*nument-1) = dR(2,:)';
                B(3,3:3:3*nument) = dR(3,:)';
                
                B(4,1:3:3*nument-2) = dR(2,:)';
                B(4,2:3:3*nument-1) = dR(1,:)';
                
                B(5,2:3:3*nument-1) = dR(3,:)';
                B(5,3:3:3*nument) = dR(2,:)';
                
                B(6,1:3:3*nument-2) = dR(3,:)';
                B(6,3:3:3*nument) = dR(1,:)';
                
                stressvect = Cmat*B*sol0(tscrtx);
                
                %set-up linear system for solving for the recovered stress
                M_int(row_index,nodes_loc) = N';
                rhs_int(row_index,:) = stressvect';
                
            end
        end
    end
    rec_sol = M_int\rhs_int;
    
    
    %calculate the error in energy norm for each element of the octuple
    for e=1:8 %for each element in the octuple
        
        e_glob = octupleList(j,e);
        
        nument = size(PHTelem(e_glob).C,1); %number of basis functions with support on current knotspan
        scrtx = PHTelem(e_glob).nodesGlobal(1:nument);
        tscrtx = reshape([3*scrtx-2; 3*scrtx-1; 3*scrtx],1,3*nument);
        nodes_loc = elementNodeLoc(e,:);
        
        Ce = PHTelem(e_glob).C;
        Ce_loc = C_loc(:,:,e);
        
        xmin = PHTelem(e_glob).vertex(1);
        xmax = PHTelem(e_glob).vertex(4);
        ymin = PHTelem(e_glob).vertex(2);
        ymax = PHTelem(e_glob).vertex(5);
        zmin = PHTelem(e_glob).vertex(3);
        zmax = PHTelem(e_glob).vertex(6);
        
        detJacobian = (xmax - xmin)*(ymax - ymin)*(zmax-zmin)/8;
        
        h1normerr = 0;
        for kk=1:numGaussW
            for jj=1:numGaussV
                for ii=1:numGaussU
                    u_hat = gaussLocationsU(ii); %gauss points on [-1, 1]
                    v_hat = gaussLocationsV(jj);
                    w_hat = gaussLocationsW(kk);
                    
                    [~, dN] = phtBasis(u_hat, v_hat, w_hat, Ce, p, q, r);
                    [~, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax);
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dN(1,:) = dN(1,:)*2/(xmax-xmin);
                    dN(2,:) = dN(2,:)*2/(ymax-ymin);
                    dN(3,:) = dN(3,:)*2/(zmax-zmin);
                    
                    % Solve for first derivatives in global coordinates
                    dN = dxdxi\dN;
                    J = abs(det(dxdxi));
                    shape = phtBasis(u_hat, v_hat, w_hat, Ce_loc, p_loc, q_loc, r_loc);
                    recovStress = shape*rec_sol(nodes_loc,:);
                    
                    %initialize the B-matrix for 3D
                    B = zeros(6,3*nument);
                    B(1,1:3:3*nument-2) = dN(1,:)';
                    B(2,2:3:3*nument-1) = dN(2,:)';
                    B(3,3:3:3*nument) = dN(3,:)';
                    
                    B(4,1:3:3*nument-2) = dN(2,:)';
                    B(4,2:3:3*nument-1) = dN(1,:)';
                    
                    B(5,2:3:3*nument-1) = dN(3,:)';
                    B(5,3:3:3*nument) = dN(2,:)';
                    
                    B(6,1:3:3*nument-2) = dN(3,:)';
                    B(6,3:3:3*nument) = dN(1,:)';
                    
                    compStress = Cmat*B*sol0(tscrtx);
                    recovStress = recovStress';
                    h1norm_glob = h1norm_glob + recovStress'*invC*recovStress*gaussWeightsU(ii)*gaussWeightsV(jj)*gaussWeightsW(kk)*detJacobian*J;
                    h1normerr = h1normerr + (recovStress-compStress)'*invC*(recovStress-compStress)*gaussWeightsU(ii)*gaussWeightsV(jj)*gaussWeightsW(kk)*detJacobian*J;
                    
                end
            end
        end
        err_glob = err_glob + h1normerr;
        estError(e_glob) = h1normerr;
    end
end

estErrorGlob = sqrt(err_glob/h1norm_glob);

local_target = targetScale/size(octupleList,1)*target_rel_error^2*h1norm_glob;
for j=1:size(octupleList,1)
    sumError = sum(estError(octupleList(j,:)));
    if sumError > local_target
        %mark the patch for refinement
        octupleRef(1,j) = sumError;
    end
end
