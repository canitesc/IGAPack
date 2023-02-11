function [nrb_new] = geometry_deg_elevate(knotU, knotV, coefs, p, q)
% elevates the geometry given in terms of knot vectors and control points
% in the NURBS toolbox format to the degrees specified by p and q
nrb = nrbmak(coefs,{knotU,knotV});
p_init = nrb.order(1)-1;
q_init = nrb.order(2)-1;
%increase polynomial order
nrb_new = nrbdegelev(nrb,[p-p_init,q-q_init]);

end

