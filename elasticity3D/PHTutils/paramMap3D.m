function [ coord, dxdxi] = paramMap3D( GIFTmesh, u_hat, v_hat, w_hat, xmin, ymin, zmin, xmax, ymax, zmax)
%maps (u_hat, v_hat, w_hat) from the reference coordinate, to physical space using the geometry
%information stored in GIFTmesh

%map (u_hat, v_hat) on [-1, 1]x[-1,1] to (xi, eta) on [xmin, xmax]x[ymin,
%ymax]

xi = u_hat*(xmax-xmin)/2+(xmax+xmin)/2;
eta = v_hat*(ymax-ymin)/2+(ymax+ymin)/2;
zeta = w_hat*(zmax-zmin)/2+(zmax+zmin)/2;

%locate the vertices of the element on the coarse mesh

for inE = 1:GIFTmesh.numberElements
    if (xi<=GIFTmesh.elementVertex(inE,4)) && (xi>=GIFTmesh.elementVertex(inE,1)) && (eta<=GIFTmesh.elementVertex(inE,5)) && (eta>=GIFTmesh.elementVertex(inE,2)) && (zeta<=GIFTmesh.elementVertex(inE,6)) && (zeta>=GIFTmesh.elementVertex(inE,3))
        e = inE;
        gxmin = GIFTmesh.elementVertex(e,1);
        gxmax = GIFTmesh.elementVertex(e,4);
        gymin = GIFTmesh.elementVertex(e,2);
        gymax = GIFTmesh.elementVertex(e,5);
        gzmin = GIFTmesh.elementVertex(e,3);
        gzmax = GIFTmesh.elementVertex(e,6);
    end
end

%map (xi, eta, zeta) on [gxmin, gxmax]x[gymin, gymax]x[gzmin, gzmax] to (uu_hat, vv_hat, ww_hat) on [-1,
%1]x[-1,1]x[-1,1];

uu_hat = (2*xi - gxmin - gxmax)/(gxmax-gxmin);
vv_hat = (2*eta - gymin - gymax)/(gymax-gymin);
ww_hat = (2*zeta - gzmin - gzmax)/(gzmax-gzmin);

nodes = GIFTmesh.elementNode(e,:);
wgts = GIFTmesh.c_net(nodes, 4);
[N, dN] = nurbshape3d_gift(uu_hat, vv_hat, ww_hat, wgts, GIFTmesh.C(:,:,e), GIFTmesh.p, GIFTmesh.q, GIFTmesh.r);
cpts = GIFTmesh.c_net(nodes,1:3);

%multiply by the jacobian of the transformation from reference
%space to the parameter space
dN(1,:) = dN(1,:)*2/(gxmax-gxmin);
dN(2,:) = dN(2,:)*2/(gymax-gymin);
dN(3,:) = dN(3,:)*2/(gzmax-gzmin);

%calculate the coordinates in the physical space
coord = N*cpts;          

%calculate the Jacobian of the transformation
dxdxi = dN*cpts;




