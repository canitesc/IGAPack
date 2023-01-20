function [ PHTmesh ] = init1DGeometryPHT(object_type, L, p, numberElements)
%creates a 1d PHT mesh of object type, polynomial degree p, numberElements elements
%includes pre-processing

PHTmesh = struct;
if strcmp(object_type, 'beam')
    rep_knot = linspace(0,1,numberElements+1);
    rep_knot = rep_knot(2:end-1);
    rep_knot = repmat(rep_knot,1,p-2);
    PHTmesh.knotVector = [zeros(1,p), linspace(0, 1, numberElements+1), ones(1,p)];
    PHTmesh.knotVector = sort([PHTmesh.knotVector, rep_knot]);
    
    lenu = length(PHTmesh.knotVector);
    numb = lenu-p-1; %number of basis functions

    %add ControlPoints at Greville abscissae
    PHTmesh.b_net = ones(numb, 2);
    for i=1:numb % for each node in the x direction
         coordx = (sum(PHTmesh.knotVector(i+1:i+p))./p);    
         PHTmesh.b_net(i,1) = L.*coordx;                   
    end
    
    PHTmesh.numberElements = numberElements;
    PHTmesh.p = p;        
    PHTmesh.L = L;
    
end

