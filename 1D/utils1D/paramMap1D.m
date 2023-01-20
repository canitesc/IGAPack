function [ coord, dxdxi] = paramMap1D( GeoMesh, u_hat, umin, umax )
%maps (u_hat, v_hat) from the reference coordinate, to physical space using the geometry
%information stored in GIFTmesh

%map (u_hat) on [-1, 1] to (xi) on [umin, umax]
xi = u_hat*(umax-umin)/2+(umax+umin)/2;

%locate the vertices of the element on the coarse mesh
for inE = 1:GeoMesh.numberElements
    if (xi<=GeoMesh.elementVertex(inE,2)) && (xi>=GeoMesh.elementVertex(inE,1)) 
        e = inE;
        gxmin = GeoMesh.elementVertex(e,1);
        gxmax = GeoMesh.elementVertex(e,2);        
    end
end

%map (xi) on [gxmin, gxmax] to (uu_hat) on [-1, 1];
uu_hat = (2*xi - gxmin - gxmax)/(gxmax-gxmin);

%evaluate the number shape function for the geometry
nodes = GeoMesh.elementNode(e,:);
wgts = GeoMesh.c_net(nodes, 2);
[N, dN] = nrbShape1DGeo(uu_hat, wgts, GeoMesh.C(:,:,e), GeoMesh.p);
cpts = GeoMesh.c_net(nodes,1);

%multiply by the jacobian of the transformation from reference
%space to the parameter space
dN = dN*2/(gxmax-gxmin);

%calculate the coordinates in the physical space
coord = N*cpts;          

%calculate the Jacobian of the transformation
dxdxi = dN*cpts;




