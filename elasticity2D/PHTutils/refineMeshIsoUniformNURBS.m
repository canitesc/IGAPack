function [nrb_new] = refineMeshIsoUniformNURBS(nrb)
knotU = nrb.knots{1};
knotV = nrb.knots{2};
p = nrb.order(1) - 1;
q = nrb.order(2) - 1;
[~, ~, new_knots_U] = kntrefine(knotU, 1, p, p-1);
[~, ~, new_knots_V] = kntrefine(knotV, 1, q, q-1);
nrb_new = nrbkntins(nrb, {new_knots_U, new_knots_V});

