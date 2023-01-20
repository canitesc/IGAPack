function [tupleRef, estErrorGlob] = recoveryDoerfler1D(PHTelem, GeoMesh, sol0, tupleList, p, targetScale)
%recover by refining the elements with largest contribution to error, not
%relative contribution
%use Bezier extraction

estError = zeros(1, length(PHTelem)); %return an array containing the estimated error in each patch and a flag on whether the element should be refined
tupleRef = zeros(1, size(tupleList,1)); %flag which patches to be refined

numGauss = p+1;
[gaussWeights, gaussLocations] = quadrature(numGauss, 'GAUSS', 1);

err_glob = 0;
h1norm_glob = 0;
toleq = 1e-10;
up = returnSuper(p);

numTuples = 0;

for j=1:size(tupleList,1)
    
    %map the superconvergent points from [-1,1] to [0..1].
    if floor(p/2)==(p/2)
        sConvPts = ([up, up+2]+1)/4;
    else
        sConvPts = ([up, up(2:end)+2]+1)/4;
    end
    
    if p==3
        p_loc = p;
        knotU = [zeros(1,p_loc+1),1/2,ones(1,p_loc+1)];
    elseif p==4
        p_loc = 5;
        knotU = [zeros(1,p_loc+1),1/2,1/2,ones(1,p_loc+1)];
    elseif p==5
        p_loc = 7;
        knotU = [zeros(1,p_loc+1),1/2,ones(1,p_loc+1)];
    elseif p==6
        p_loc = 6;
        knotU = [zeros(1,p_loc+1),1/2,ones(1,p_loc+1)];
    elseif p==7
        p_loc = 7;
        knotU = [zeros(1,p_loc+1),1/2,ones(1,p_loc+1)];
    end
    
    
    lenU = length(knotU)-p_loc-1;  %number of basis functions in the u direction
    
    elementCounterLoc = 0;
    elementVertexLoc = zeros(2, 2);
    elementNodeLoc = zeros(2, p_loc+1);
    
    for ii=1:length(knotU)-1
        if (abs(knotU(ii+1)-knotU(ii))>toleq)
            elementCounterLoc = elementCounterLoc + 1;
            elementVertexLoc(elementCounterLoc, :) = [knotU(ii), knotU(ii+1)];
            
            currow = ii-p_loc:ii;
            elementNodeLoc(elementCounterLoc,:)=currow;
        end
    end
    
    M_int = sparse(length(sConvPts),lenU); %the interpolation matrix
    rhs_int = zeros(length(sConvPts),1); %interpolation RHS
    C_loc = bezierExtraction(knotU, p_loc);
    
    for i=1:length(sConvPts)
        ucoord = sConvPts(i);
        if (ucoord >= 0) && (ucoord <= 0.5)
            e = 1;
            e_glob = tupleList(j,1);
            u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
        else
            e = 2;
            e_glob = tupleList(j,2);
            u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
        end
        nodes_loc = elementNodeLoc(e,:);
        
        umin = PHTelem(e_glob).vertex(1);
        umax = PHTelem(e_glob).vertex(2);
        
        N = phtBasis1D(u_hat, C_loc(:,:,e), p_loc);
        Ce = PHTelem(e_glob).C;
        
        [ ~, dxdxi] = paramMap1D( GeoMesh, u_hat, umin, umax );
        [~, dR] = phtBasis1D( u_hat, Ce, p );
        
        nument = size(PHTelem(e_glob).C,1); %number of basis functions with support on current knotspan
        sctrx = PHTelem(e_glob).nodes(1:nument);
        
        %multiply by the jacobian of the transformation from reference
        %space to the parameter space
        dR = dR*2/(umax-umin);
        dR = dxdxi\dR;
        
        M_int(i,nodes_loc) = N;
        
        rhs_int(i) = dR*sol0(sctrx);
    end
    
    rec_sol = M_int\rhs_int;
    
    %calculate the error in H1 norm for each element of the patch
    for e=1:2 %for each element in the patch
        
        e_glob = tupleList(j,e);
        nument = size(PHTelem(e_glob).C,1); %number of basis functions with support on current knotspan
        sctrx = PHTelem(e_glob).nodes(1:nument);
        nodes_loc = elementNodeLoc(e,:);
        
        Ce = PHTelem(e_glob).C;
        Ce_loc = C_loc(:,:,e);
        
        umin = PHTelem(e_glob).vertex(1);
        umax = PHTelem(e_glob).vertex(2);
        
        detJacobian = (umax-umin)/2;
        
        h1normerr = 0;
        for q=1:numGauss
            u_hat = gaussLocations(q); %gauss points on [-1, 1]
            
            [~, dN] = phtBasis1D(u_hat, Ce, p);
            [~, dxdxi] = paramMap1D(GeoMesh, u_hat, umin, umax);
            
            %multiply by the jacobian of the transformation from reference
            %space to the parameter space
            dN = dN*2/(umax-umin);
            
            % Solve for first derivatives in global coordinates
            dN = dxdxi\dN;
            
            shape = phtBasis1D(u_hat, Ce_loc, p_loc);
            
            recovDeriv = shape*rec_sol(nodes_loc);
            compDeriv = dN*sol0(sctrx);
            
            h1norm_glob = h1norm_glob + recovDeriv^2*gaussWeights(q)*detJacobian;
            h1normerr = h1normerr + (recovDeriv-compDeriv)^2*gaussWeights(q)*detJacobian;
            
        end
        
        err_glob = err_glob + h1normerr;
        estError(1,e_glob) = h1normerr;
    end
   
end

numTuples = size(tupleList,1);
estErrorGlob = sqrt(err_glob/h1norm_glob);


tupleCounter = 0;
tupleErrors = zeros(1, numTuples);
tupleIndex = zeros(numTuples,2);
for j=1:size(tupleList,1)
    tupleCounter = tupleCounter + 1;
    tupleErrors(tupleCounter) = sum(estError(tupleList(j,:)));
    tupleIndex(tupleCounter,:) = [1, j];
end

[tupleErrors, index] = sort(tupleErrors,'descend');
tupleIndex = tupleIndex(index, :);

targetError = targetScale*sqrt(err_glob);

sumError = 0;
fudge = 1e-10; %parameter added to ensure uniform refinement when targetScale = 1
%note there may be some elements where the estimated error is less than eps
for i=1:numTuples
    sumError = sumError + tupleErrors(i);
    %i
    %tupleIndex(i,2)
    tupleRef(tupleIndex(i,2)) = 1;
    if sqrt(sumError) > targetError+fudge
        disp(['Refining ', num2str(i), ' quads with largest error out of ', num2str(numTuples)])
        break
    end
end
