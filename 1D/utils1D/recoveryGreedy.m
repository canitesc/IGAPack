function [estError,patchRef,estErrorGlob]=recoveryGreedy(PHTmesh,elementNodes,displacements,target_rel_error,patchList,C)
%recover by refining the elements with largest contribution to error, not
%relative contribution
%use Bezier extraction

estError = zeros(2, PHTmesh.numberElements); %return an array containing the estimated error in each patch and a flag on whether the element should be refined
patchRef = zeros(1, size(patchList,1)); %flag which patches to be refined
knotVector = PHTmesh.knotVector;
b_net = PHTmesh.b_net;
p = PHTmesh.p;
numGauss = p+1;
[gaussWeights, gaussLocations] = quadrature(numGauss, 'GAUSS', 1);

err_glob = 0;
h1norm_glob = 0;
toleq = 1e-10;
%figure

for j=1:size(patchList,1)
    up = returnSuper(p);
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
    
    %add ControlPoints at Greville abscissae
    b_net_loc = ones(lenU, 2);
    for i=1:lenU % for each node in the x direction
        coordx = (sum(knotU(i+1:i+p_loc))./p_loc);
        b_net_loc(i,1) = coordx;
    end
    
    
    for ii=1:length(knotU)-1
        if (abs(knotU(ii+1)-knotU(ii))>toleq)
            elementCounterLoc = elementCounterLoc + 1;
            elementVertexLoc(elementCounterLoc, :) = [knotU(ii), knotU(ii+1)];
            tcount = 0;
            currow = zeros(1, (p_loc+1));
            %now we add the nodes from i-p...i in the u direction
            for t1 = ii-p_loc:ii
                tcount = tcount + 1;
                currow(tcount) = t1;
            end
            
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
            e_glob = patchList(j,1);
            u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
        else
            e = 2;
            e_glob = patchList(j,2);
            u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
        end
        nodes = elementNodes(e_glob,:);
        nodes_loc = elementNodeLoc(e,:);
        
        Ce_loc = C_loc(:,:,e);
        
        cpts_loc = b_net_loc(nodes_loc, 1);
        w_loc = b_net_loc(nodes_loc, 2);
        [N] = basisFun(u_hat, cpts_loc, w_loc, Ce_loc, p_loc);
        
        Ce = C(:,:,e_glob);
        cpts = b_net(nodes, 1);
        w = b_net(nodes, 2);
        [~, ~, ~, dN] = basisFun(u_hat, cpts, w, Ce, p);
        
        M_int(i,nodes_loc) = N';
        rhs_int(i) = dN*displacements(nodes);
    end
    
    rec_sol = M_int\rhs_int;
    
%     %plot the estimated error
%     numPoints = 50;
%     plotPoints = linspace(0,1,numPoints);
%     
%     coordPoints = zeros(1,numPoints);
%     errorPoints = zeros(1,numPoints);
%     
%     %shapePoints = zeros(1,numPlot);
%     for i=1:numPoints
%         ucoord = plotPoints(i);
%         if (ucoord >= 0) && (ucoord <= 0.5)
%             e = 1;
%             e_glob = patchList(j,1);
%             u_hat = 4*ucoord-1; %map [0,0.5]-> [-1,1]
%         else
%             e = 2;
%             e_glob = patchList(j,2);
%             u_hat = 4*ucoord-3; %map [0.5,1] -> [-1,1]
%         end
%         
%         Ce = C(:,:,e_glob);
%         Ce_loc = C_loc(:,:,e);
%         
%         nodes_glob = elementNodes(e_glob,:);
%         cpts = b_net(nodes_glob, 1);
%         w = b_net(nodes_glob, 2);
%         
%         nodes_loc = elementNodeLoc(e,:);
%         cpts_loc = b_net_loc(nodes_loc, 1);
%         w_loc = b_net_loc(nodes_loc, 2);
%         
%         [~,coord,~,dN] = basisFun(u_hat,cpts,w,Ce,p);
%         
%         shape = basisFun(u_hat, cpts_loc, w_loc, Ce_loc, p_loc);
%         
%         recovDeriv = shape*rec_sol(nodes_loc);
%         compDeriv = dN*displacements(nodes_glob);
%         coordPoints(i) = coord;
%         errorPoints(i) = (recovDeriv-compDeriv);
%     end
    
%     plot(coordPoints,errorPoints,'-b',coordPoints,0,'-k')
%     hold on
    
    %calculate the error in H1 norm for each element of the patch
    for e=1:2 %for each element in the patch
        if e==1
            e_glob = patchList(j,1);
        else
            e_glob = patchList(j,2);
        end
        nodes_glob = elementNodes(e_glob,:);
        nodes_loc = elementNodeLoc(e,:);
        
        cpts_loc = b_net_loc(nodes_loc, 1);
        w_loc = b_net_loc(nodes_loc, 2);
        
        cpts = b_net(nodes_glob, 1);
        w = b_net(nodes_glob, 2);
        
        Ce = C(:,:,e_glob);
        Ce_loc = C_loc(:,:,e);
        
        length_element=knotVector(nodes_glob(end)+1)-knotVector(nodes_glob(end));
        detJacobian=length_element/2;
        
        h1normerr = 0;
        for q=1:numGauss
            u_hat = gaussLocations(q); %gauss points on [-1, 1]
            
            [~,coord,~,dN] = basisFun(u_hat,cpts,w,Ce,p);
            shape = basisFun(u_hat, cpts_loc, w_loc, Ce_loc, p_loc);
            
            recovDeriv = shape*rec_sol(nodes_loc);
            compDeriv = dN*displacements(nodes_glob);
            coordPoints(i) = coord;
            
            h1norm_glob = h1norm_glob + recovDeriv^2*gaussWeights(q)*detJacobian;
            h1normerr = h1normerr + (recovDeriv-compDeriv)^2*gaussWeights(q)*detJacobian;
            
        end
        
        err_glob = err_glob + h1normerr;
        estError(1,e_glob) = h1normerr;
    end
end

estErrorGlob = sqrt(err_glob/h1norm_glob);
local_target = 1/PHTmesh.numberElements*target_rel_error^2*h1norm_glob;
for j=1:size(patchList,1)
    i = patchList(j,1); %pick the first element
    
    if estError(1,i) > local_target %mark the elements with error above threshold for refinement
        nodes_glob = elementNodes(i,:);
        estError(2,i) = 1;
        patchRef(j) = 1;
        %line([knotVector(nodes_glob(end)+1),knotVector(nodes_glob(end))],[0,0],'LineWidth',3,'Color','red');
    end
    i = patchList(j,2); %pick the second element
    if estError(1,i) > local_target %mark the elements with error above threshold for refinement
        nodes_glob = elementNodes(i,:);
        estError(2,i) = 1;
        patchRef(j) = 1;
        %line([knotVector(nodes_glob(end)+1),knotVector(nodes_glob(end))],[0,0],'LineWidth',3,'Color','red');
    end
end
% coordPoints = linspace(0,1,1001);
% plot(coordPoints, target_rel_error*sqrt(h1norm_glob), '-r', coordPoints, -target_rel_error*sqrt(h1norm_glob), '-r')
% title('Estimated Error')
% hold off
% drawnow
