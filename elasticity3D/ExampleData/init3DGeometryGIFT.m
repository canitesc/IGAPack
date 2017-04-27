function [ GIFTmesh] = init3DGeometryGIFT(object_type,L,W,H)
%creates a 2d GIFT mesh and associated control points, knotvectors, Bezier
%extraction operators, etc.

if strcmp(object_type, 'cube')
    
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    %initialize geometry on the coarsest mesh
    coefs(1:3,1,1,1) = [0; 0; 0];
    coefs(1:3,1,2,1) = [0; W; 0];
    coefs(1:3,2,1,1) = [L; 0; 0];
    coefs(1:3,2,2,1) = [L; W; 0];
    coefs(1:3,1,1,2) = [0; 0; H];
    coefs(1:3,1,2,2) = [0; W; H];
    coefs(1:3,2,1,2) = [L; 0; H];
    coefs(1:3,2,2,2) = [L; W; H];
    coefs(4,1,1,1) = 1;
    coefs(4,1,2,1) = 1;
    coefs(4,2,1,1) = 1;
    coefs(4,2,2,1) = 1;
    coefs(4,1,1,2) = 1;
    coefs(4,1,2,2) = 1;
    coefs(4,2,1,2) = 1;
    coefs(4,2,2,2) = 1;
    
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    knotW = [0 0 1 1];
    
    [ GIFTmesh ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    
elseif strcmp(object_type, 'horseshoe')
    %use the C1 horseshoe model obtained by extruding/revolving the plate with a hole
    numberElementsU = 5;
    numberElementsV = 2;
    numberElementsW = 1;
    p = 2;
    q = 2;
    r = 2;
    
    %knotU = [0 0 0 0.5000 1 1 1];
    %knotV = [0 0 0 1 1 1];
    %knotW = [0, 0, 0, 1/6, 1/3, 4/9, 5/9, 2/3, 5/6, 1, 1, 1];
    
    %load the coefficients from the variable file
    load('horse.mat')
    
    [ GIFTmesh ] = genGIFTmesh3D(horse.knots{1}, horse.knots{2}, horse.knots{3}, horse.coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
elseif strcmp(object_type, 'blade')
    numberElementsU = 10;
    numberElementsV = 20;
    numberElementsW = 1;
    p = 2;
    q = 2;
    r = 1;
    load('blade.mat')
    [ GIFTmesh ] = genGIFTmesh3D(blade.knots{1}, blade.knots{2}, blade.knots{3}, blade.coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);

    
elseif strcmp(object_type, 'hemisphere')
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p=2;
    q=2;
    r=1;
    load('hemisphere.mat')
    [ GIFTmesh ] = genGIFTmesh3D(hemisphere.knots{1}, hemisphere.knots{2}, hemisphere.knots{3}, hemisphere.coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    
elseif strcmp(object_type, 'cylinder')
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p=1;
    q=2;
    r=1;
    cylinder = [];
    load('cylinder.mat')
    %     p=2;
    %     r=2;
    %     cylinder = nrbdegelev(cylinder, [p,q,r]-(cylinder.order-1));
    
    [ GIFTmesh ] = genGIFTmesh3D(cylinder.knots{1}, cylinder.knots{2}, cylinder.knots{3}, cylinder.coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
elseif strcmp(object_type, 'cube_with_hole')
    numberElementsU = 2;
    numberElementsV = 2;
    numberElementsW = 1;
    p=2;
    q=2;
    r=2;
    cube_with_hole = [];
    load('cube_with_hole.mat')
   % figure;
   % nrbkntplot(cube_with_hole)
    [ GIFTmesh ] = genGIFTmesh3D(cube_with_hole.knots{1}, cube_with_hole.knots{2}, cube_with_hole.knots{3}, cube_with_hole.coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
end




