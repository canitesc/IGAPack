%function [stiffness,force] = stiffness1d(numuniqelem, elementNodes, numb, b_net, knotVector, p)
function [stiffness,force] = stiffness1d(PHTmesh, elementNodes, C)
% calculate stiffness matrix for 1d problem


p = PHTmesh.p;
knotVector = PHTmesh.knotVector;
numBasis = length(knotVector)-p-1; %number of basis functions
stiffness=sparse(numBasis,numBasis);
force=zeros(numBasis,1);
b_net = PHTmesh.b_net;
numberElements = PHTmesh.numberElements;

% stiffness matrix
[gaussWeights, gaussLocations] = quadrature(p+1, 'GAUSS', 1);

for k=1:numberElements
    elementDof=elementNodes(k,:);   
    ndof=p+1;
    length_element=knotVector(elementDof(end)+1)-knotVector(elementDof(end));
    detJacobian=length_element/2;
    Ce = C(:,:,k);   
    cpts = PHTmesh.b_net(elementDof,1);
    w = PHTmesh.b_net(elementDof, 2);
    %invJacobian=1/detJacobian;
    for q=1:size(gaussWeights,1)
        pt=gaussLocations(q,:);
        [N, coord, ~, dN, detJ] = basisFun(pt, cpts, w, Ce, p);               
        [a0, a1] = coefPDE(coord);
        
        % B matrix
        B=zeros(1, ndof);
        M=zeros(1, ndof);
        M(1:ndof) = N'; 
        B(1:ndof) = dN';
        
        % K       
        stiffness(elementDof,elementDof)=...
            stiffness(elementDof,elementDof)+...
            (B'*B+a1*M'*B+a0*M'*M)*gaussWeights(q)*detJacobian*detJ;
        [~, ~, P] = exact_sol1d(coord);
        force(elementDof)=force(elementDof)+N'*P*detJacobian*gaussWeights(q)*detJ;
       
    end
end