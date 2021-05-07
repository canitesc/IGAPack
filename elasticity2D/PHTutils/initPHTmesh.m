function [ PHTelem, dimBasis ] = initPHTmesh( p,q )
%initialize the PHT geometry on coarse mesh

%start with num_elements_u, num_elements_v uniformly distributed knot spans
%in each parametric direction
num_elements_u = 1;
num_elements_v = 1;
dimBasis = (p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1)

knotU = [zeros(1,p+1), (1:num_elements_u)/num_elements_u, ones(1,p)];
knotV = [zeros(1,q+1), (1:num_elements_v)/num_elements_v, ones(1,q)];

%repeat the interior knots p-1 times
rep_knotU = linspace(0,1,num_elements_u+1);
rep_knotU = rep_knotU(2:end-1);
rep_knotU = repmat(rep_knotU,1,p-2);

rep_knotV = linspace(0,1,num_elements_v+1);
rep_knotV = rep_knotV(2:end-1);
rep_knotV = repmat(rep_knotV,1,q-2);

knotU = sort([knotU, rep_knotU]);
knotV = sort([knotV, rep_knotV]);


[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);

PHTelem = struct;

%define the 1st element
PHTelem.parent = [];
PHTelem.children = [];
PHTelem.C = kron(C_v(:,:,1),C_u(:,:,1));
PHTelem.vertex = [0, 0, 1, 1];
PHTelem.nodes = 1:(p+1)*(q+1);
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];
PHTelem.level = 0;

[ PHTelem, dimBasis ] = crossInsert( PHTelem, 1, dimBasis, p, q );
%[ PHTelem, dimBasis ] = crossInsert( PHTelem, 5, dimBasis, p, q );
