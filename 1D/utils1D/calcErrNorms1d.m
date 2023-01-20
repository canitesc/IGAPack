function [l2normerr,h1normerr]=calcErrNorms1d(PHTmesh,elementNodes,displacements,C)
% calculate the L2 norm of the error for NURBS basis functions
% GDof: total number of degrees of freedom of
% the problem
p = PHTmesh.p;
numGauss = p+1;
[gaussWeights, gaussLocations] = quadrature(numGauss, 'GAUSS', 1);

knotVector = PHTmesh.knotVector;
b_net = PHTmesh.b_net;

l2norm = 0;
h1norm = 0;
l2normerr = 0;
h1normerr = 0;

%loop over the elements
for k=1:PHTmesh.numberElements
    
    nodes = elementNodes(k,:);  
    cpts = b_net(nodes,1);
    w = b_net(nodes,2);
    
    length_element=knotVector(nodes(end)+1)-knotVector(nodes(end));
    detJacobian=length_element/2;
    Ce = C(:,:,k);    
    
    for q=1:numGauss
        pt = gaussLocations(q);
        
        [shape,coord,~,shgradg] = basisFun(pt,cpts,w,Ce,p);
        compSol = shape*displacements(nodes);
        compDeriv = shgradg*displacements(nodes);
        
        [analSol, analDeriv] = exact_sol1d(coord);
        
        l2norm = l2norm + analSol^2*gaussWeights(q)*detJacobian;
        h1norm = h1norm + analDeriv^2*gaussWeights(q)*detJacobian;
        
        l2normerr = l2normerr + (analSol-compSol)^2*gaussWeights(q)*detJacobian;    
        h1normerr = h1normerr + (analDeriv-compDeriv)^2*gaussWeights(q)*detJacobian;        
    end
end

l2normerr = sqrt(l2normerr)/sqrt(l2norm);
h1normerr = sqrt(h1normerr)/sqrt(h1norm);
