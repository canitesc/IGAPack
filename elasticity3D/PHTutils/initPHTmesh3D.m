function [ PHTelem, dimBasis ] = initPHTmesh3D( p,q,r )
%initialize the PHT geometry on coarse mesh

%start with num_elements_u, num_elements_v, num_elements_w uniformly distributed knot spans
%in each parametric direction
num_elements_u = 1;
num_elements_v = 1;
num_elements_w = 1;
dimBasis = (p+1)*(q+1)*(r+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1)*(r+1)

knotU = [zeros(1,p+1), (1:num_elements_u)/num_elements_u, ones(1,p)];
knotV = [zeros(1,q+1), (1:num_elements_v)/num_elements_v, ones(1,q)];
knotW = [zeros(1,r+1), (1:num_elements_w)/num_elements_w, ones(1,r)];

%repeat the interior knots p-1 times
rep_knotU = linspace(0,1,num_elements_u+1);
rep_knotU = rep_knotU(2:end-1);
rep_knotU = repmat(rep_knotU,1,p-2);

rep_knotV = linspace(0,1,num_elements_v+1);
rep_knotV = rep_knotV(2:end-1);
rep_knotV = repmat(rep_knotV,1,q-2);

rep_knotW = linspace(0,1,num_elements_w+1);
rep_knotW = rep_knotW(2:end-1);
rep_knotW = repmat(rep_knotW,1,r-2);

knotU = sort([knotU, rep_knotU]);
knotV = sort([knotV, rep_knotV]);
knotW = sort([knotW, rep_knotW]);


[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);
[C_w, ~] = bezierExtraction(knotW,r);

PHTelem = struct;

%define the 1st element
PHTelem.parent = [];
PHTelem.children = [];
PHTelem.C = kron(kron(C_w(:,:,1),C_v(:,:,1)),C_u(:,:,1));
PHTelem.vertex = [0, 0, 0, 1, 1, 1];
PHTelem.nodes = 1:(p+1)*(q+1)*(r+1);
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];
PHTelem.neighbor_front = [];
PHTelem.neighbor_back = [];
PHTelem.neighbor_up_left = [];
PHTelem.neighbor_down_left = [];
PHTelem.neighbor_up_right = [];
PHTelem.neighbor_down_right = [];
PHTelem.neighbor_up_front = [];
PHTelem.neighbor_down_front = [];
PHTelem.neighbor_up_back = [];
PHTelem.neighbor_down_back = [];
PHTelem.neighbor_left_front = [];
PHTelem.neighbor_right_front = [];
PHTelem.neighbor_left_back = [];
PHTelem.neighbor_right_back = [];

PHTelem.level = 0;

[ PHTelem, dimBasis ] = crossInsert3D( PHTelem, 1, dimBasis, p, q ,r );

