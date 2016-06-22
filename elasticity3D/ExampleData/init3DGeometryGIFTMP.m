function [ GIFTmesh] = init3DGeometryGIFTMP(object_type,L,W,H,numPatches,numPatchesU,numPatchesV)
%creates a 2d GIFT mesh and associated control points, knotvectors, Bezier
%extraction operators, etc.

if strcmp(object_type, 'cube')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,numPatches+1);
    
    for patchIndex = 1:numPatches
        
        %set the dimensions of the patch
        patchMinX = xVertices(patchIndex);
        patchMaxX = xVertices(patchIndex+1);
        patchMinY = 0;
        patchMaxY = W;
        patchMinZ = 0;
        patchMaxZ = H;
        
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1,1) = [patchMinX; patchMinY; patchMinZ];
        coefs(1:3,1,2,1) = [patchMinX; patchMaxY; patchMinZ];
        coefs(1:3,2,1,1) = [patchMaxX; patchMinY; patchMinZ];
        coefs(1:3,2,2,1) = [patchMaxX; patchMaxY; patchMinZ];
        coefs(1:3,1,1,2) = [patchMinX; patchMinY; patchMaxZ];
        coefs(1:3,1,2,2) = [patchMinX; patchMaxY; patchMaxZ];
        coefs(1:3,2,1,2) = [patchMaxX; patchMinY; patchMaxZ];
        coefs(1:3,2,2,2) = [patchMaxX; patchMaxY; patchMaxZ];
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
        
        [ GIFTmesh{patchIndex} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    end
elseif strcmp(object_type, 'cubeUV')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    %divide the patches along the x and y directions
    xVertices = linspace(0,L,numPatchesU+1);
    yVertices = linspace(0,W,numPatchesV+1);
    patchIndex = 0;
    for patchIndexV = 1:numPatchesV
        for patchIndexU = 1:numPatchesU
            patchIndex = patchIndex+1;
            %set the dimensions of the patch
            patchMinX = xVertices(patchIndexU);
            patchMaxX = xVertices(patchIndexU+1);
            patchMinY = yVertices(patchIndexV);
            patchMaxY = yVertices(patchIndexV+1);
            patchMinZ = 0;
            patchMaxZ = H;
            
            %initialize geometry on coarsest mesh
            coefs(1:3,1,1,1) = [patchMinX; patchMinY; patchMinZ];
            coefs(1:3,1,2,1) = [patchMinX; patchMaxY; patchMinZ];
            coefs(1:3,2,1,1) = [patchMaxX; patchMinY; patchMinZ];
            coefs(1:3,2,2,1) = [patchMaxX; patchMaxY; patchMinZ];
            coefs(1:3,1,1,2) = [patchMinX; patchMinY; patchMaxZ];
            coefs(1:3,1,2,2) = [patchMinX; patchMaxY; patchMaxZ];
            coefs(1:3,2,1,2) = [patchMaxX; patchMinY; patchMaxZ];
            coefs(1:3,2,2,2) = [patchMaxX; patchMaxY; patchMaxZ];
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
            
            [ GIFTmesh{patchIndex} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
        end
    end
elseif strcmp(object_type, 'Tbeam')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    Lmax = 10;
    Lmid = 5;
    w1 = 0.6;
    h1 = 0.6;
    w2 = 1;
    h2 = 1;
    
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    knotW = [0 0 1 1];
    
    %initialize geometry on the each patch
    %set the dimensions of the patch
    patchMinX = [0, 0, 0, 0, Lmid];
    patchMaxX = [Lmid, Lmid, Lmid, Lmid, Lmax];
    patchMinY = [0, w1, 0, w1, 0];
    patchMaxY = [w1, w2, w1, w2, w1];
    patchMinZ = [0, 0, h1, h1, 0];
    patchMaxZ = [h1, h1, h2, h2, h1];
    
    for indexPatch=1:5
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1,1) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,2,1) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,1,1) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,2,1) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,1,2) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,1,2,2) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,1,2) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,2,2) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(4,1,1,1) = 1;
        coefs(4,1,2,1) = 1;
        coefs(4,2,1,1) = 1;
        coefs(4,2,2,1) = 1;
        coefs(4,1,1,2) = 1;
        coefs(4,1,2,2) = 1;
        coefs(4,2,1,2) = 1;
        coefs(4,2,2,2) = 1;
        [ GIFTmesh{indexPatch} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    end
elseif strcmp(object_type, 'ConnectingRod')
    load('ConnRod.mat')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    [ GIFTmesh{1} ] = genGIFTmesh3D(solid1A.knots{1}, solid1A.knots{2}, solid1A.knots{3}, solid1A.coefs, solid1A.order(1)-1, solid1A.order(2)-1, solid1A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{2} ] = genGIFTmesh3D(solid1B.knots{1}, solid1B.knots{2}, solid1B.knots{3}, solid1B.coefs, solid1B.order(1)-1, solid1B.order(2)-1, solid1B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{3} ] = genGIFTmesh3D(solid1C.knots{1}, solid1C.knots{2}, solid1C.knots{3}, solid1C.coefs, solid1C.order(1)-1, solid1C.order(2)-1, solid1C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{4} ] = genGIFTmesh3D(solid2A.knots{1}, solid2A.knots{2}, solid2A.knots{3}, solid2A.coefs, solid2A.order(1)-1, solid2A.order(2)-1, solid2A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{5} ] = genGIFTmesh3D(solid2B.knots{1}, solid2B.knots{2}, solid2B.knots{3}, solid2B.coefs, solid2B.order(1)-1, solid2B.order(2)-1, solid2B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{6} ] = genGIFTmesh3D(solid2C.knots{1}, solid2C.knots{2}, solid2C.knots{3}, solid2C.coefs, solid2C.order(1)-1, solid2C.order(2)-1, solid2C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{7} ] = genGIFTmesh3D(solid3A.knots{1}, solid3A.knots{2}, solid3A.knots{3}, solid3A.coefs, solid3A.order(1)-1, solid3A.order(2)-1, solid3A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{8} ] = genGIFTmesh3D(solid3B.knots{1}, solid3B.knots{2}, solid3B.knots{3}, solid3B.coefs, solid3B.order(1)-1, solid3B.order(2)-1, solid3B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{9} ] = genGIFTmesh3D(solid3C.knots{1}, solid3C.knots{2}, solid3C.knots{3}, solid3C.coefs, solid3C.order(1)-1, solid3C.order(2)-1, solid3C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{10} ] = genGIFTmesh3D(solid4.knots{1}, solid4.knots{2}, solid4.knots{3}, solid4.coefs, solid4.order(1)-1, solid4.order(2)-1, solid4.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{11} ] = genGIFTmesh3D(solidT1B.knots{1}, solidT1B.knots{2}, solidT1B.knots{3}, solidT1B.coefs, solidT1B.order(1)-1, solidT1B.order(2)-1, solidT1B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{12} ] = genGIFTmesh3D(solidT2B.knots{1}, solidT2B.knots{2}, solidT2B.knots{3}, solidT2B.coefs, solidT2B.order(1)-1, solidT2B.order(2)-1, solidT2B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{13} ] = genGIFTmesh3D(solidT3B.knots{1}, solidT3B.knots{2}, solidT3B.knots{3}, solidT3B.coefs, solidT3B.order(1)-1, solidT3B.order(2)-1, solidT3B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{14} ] = genGIFTmesh3D(solidT4B.knots{1}, solidT4B.knots{2}, solidT4B.knots{3}, solidT4B.coefs, solidT4B.order(1)-1, solidT4B.order(2)-1, solidT4B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{15} ] = genGIFTmesh3D(solidT5B.knots{1}, solidT5B.knots{2}, solidT5B.knots{3}, solidT5B.coefs, solidT5B.order(1)-1, solidT5B.order(2)-1, solidT5B.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{16} ] = genGIFTmesh3D(solidT1A.knots{1}, solidT1A.knots{2}, solidT1A.knots{3}, solidT1A.coefs, solidT1A.order(1)-1, solidT1A.order(2)-1, solidT1A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{17} ] = genGIFTmesh3D(solidT2A.knots{1}, solidT2A.knots{2}, solidT2A.knots{3}, solidT2A.coefs, solidT2A.order(1)-1, solidT2A.order(2)-1, solidT2A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{18} ] = genGIFTmesh3D(solidT3A.knots{1}, solidT3A.knots{2}, solidT3A.knots{3}, solidT3A.coefs, solidT3A.order(1)-1, solidT3A.order(2)-1, solidT3A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{19} ] = genGIFTmesh3D(solidT4A.knots{1}, solidT4A.knots{2}, solidT4A.knots{3}, solidT4A.coefs, solidT4A.order(1)-1, solidT4A.order(2)-1, solidT4A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{20} ] = genGIFTmesh3D(solidT5A.knots{1}, solidT5A.knots{2}, solidT5A.knots{3}, solidT5A.coefs, solidT5A.order(1)-1, solidT5A.order(2)-1, solidT5A.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{21} ] = genGIFTmesh3D(solidT1C.knots{1}, solidT1C.knots{2}, solidT1C.knots{3}, solidT1C.coefs, solidT1C.order(1)-1, solidT1C.order(2)-1, solidT1C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{22} ] = genGIFTmesh3D(solidT2C.knots{1}, solidT2C.knots{2}, solidT2C.knots{3}, solidT2C.coefs, solidT2C.order(1)-1, solidT2C.order(2)-1, solidT2C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{23} ] = genGIFTmesh3D(solidT3C.knots{1}, solidT3C.knots{2}, solidT3C.knots{3}, solidT3C.coefs, solidT3C.order(1)-1, solidT3C.order(2)-1, solidT3C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{24} ] = genGIFTmesh3D(solidT4C.knots{1}, solidT4C.knots{2}, solidT4C.knots{3}, solidT4C.coefs, solidT4C.order(1)-1, solidT4C.order(2)-1, solidT4C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{25} ] = genGIFTmesh3D(solidT5C.knots{1}, solidT5C.knots{2}, solidT5C.knots{3}, solidT5C.coefs, solidT5C.order(1)-1, solidT5C.order(2)-1, solidT5C.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    
    
    
elseif strcmp(object_type, 'crackcube')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    knotW = [0 0 1 1];
    
    %initialize geometry on the each patch
    %set the dimensions of the patch
    patchMinX = [-1/2, 0, 0, -1/2];
    patchMaxX = [0, 1/2, 1/2, 0];
    patchMinY = [-1/2, -1/2, -1/2, -1/2];
    patchMaxY = [1/2, 1/2, 1/2, 1/2];
    patchMinZ = [-1/2, -1/2, 0, 0];
    patchMaxZ = [0, 0, 1/2, 1/2];
    
    for indexPatch=1:4
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1,1) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,2,1) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,1,1) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,2,1) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,1,2) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,1,2,2) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,1,2) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,2,2) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(4,1,1,1) = 1;
        coefs(4,1,2,1) = 1;
        coefs(4,2,1,1) = 1;
        coefs(4,2,2,1) = 1;
        coefs(4,1,1,2) = 1;
        coefs(4,1,2,2) = 1;
        coefs(4,2,1,2) = 1;
        coefs(4,2,2,2) = 1;
        [ GIFTmesh{indexPatch} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    end
    
elseif strcmp(object_type, '2dSlitCrack')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 1;
    r = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    knotW = [0 0 1 1];
    
    %initialize geometry on the each patch
    %set the dimensions of the patch
    patchMinX = [0, 1];
    patchMaxX = [1, 2];
    patchMinY = [0, 0];
    patchMaxY = [1, 1];
    patchMinZ = [0, 0];
    patchMaxZ = [1, 1];
    
    for indexPatch=1:2
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1,1) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,2,1) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,1,1) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,2,2,1) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMinZ(indexPatch)];
        coefs(1:3,1,1,2) = [patchMinX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,1,2,2) = [patchMinX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,1,2) = [patchMaxX(indexPatch); patchMinY(indexPatch); patchMaxZ(indexPatch)];
        coefs(1:3,2,2,2) = [patchMaxX(indexPatch); patchMaxY(indexPatch); patchMaxZ(indexPatch)];
        coefs(4,1,1,1) = 1;
        coefs(4,1,2,1) = 1;
        coefs(4,2,1,1) = 1;
        coefs(4,2,2,1) = 1;
        coefs(4,1,1,2) = 1;
        coefs(4,1,2,2) = 1;
        coefs(4,2,1,2) = 1;
        coefs(4,2,2,2) = 1;
        
        
        [ GIFTmesh{indexPatch} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
    end
    
elseif strcmp(object_type, 'pennycrack')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 1;
    q = 2;
    r = 1;
    
    knotU = [0 0 1 1];
    knotV =  [0 0 0 1 1 1];
    knotW = [0, 0, 1, 1];
    rad_int = [0, 1];
    rad_ext = [1, 4];
    z_min = 0;
    z_max = H;
    
    
    %initialize geometry on coarsest meshes
    for indexPatch=1:2
        coefs(1:3,1,1,1) = [rad_int(indexPatch)*[1;0]; z_min];
        coefs(1:3,1,2,1) = [rad_int(indexPatch)*[sqrt(2)/2; sqrt(2)/2]; z_min];
        coefs(1:3,1,3,1) = [rad_int(indexPatch)*[0;1]; z_min];
        coefs(1:3,2,1,1) = [rad_ext(indexPatch)*[1;0]; z_min];
        coefs(1:3,2,2,1) = [rad_ext(indexPatch)*[sqrt(2)/2;sqrt(2)/2]; z_min];
        coefs(1:3,2,3,1) = [rad_ext(indexPatch)*[0;1]; z_min];
        coefs(1:3,1,1,2) = [rad_int(indexPatch)*[1;0]; z_max];
        coefs(1:3,1,2,2) = [rad_int(indexPatch)*[sqrt(2)/2; sqrt(2)/2]; z_max* (sqrt(2)/2)];
        coefs(1:3,1,3,2) = [rad_int(indexPatch)*[0;1]; z_max];
        coefs(1:3,2,1,2) = [rad_ext(indexPatch)*[1;0]; z_max];
        coefs(1:3,2,2,2) = [rad_ext(indexPatch)*[sqrt(2)/2;sqrt(2)/2]; z_max*(sqrt(2)/2)];
        coefs(1:3,2,3,2) = [rad_ext(indexPatch)*[0;1]; z_max];
        coefs(4,1,1,1) = 1;
        coefs(4,1,2,1) = sqrt(2)/2;
        coefs(4,1,3,1) = 1;
        coefs(4,2,1,1) = 1;
        coefs(4,2,2,1) = sqrt(2)/2;
        coefs(4,2,3,1) = 1;
        coefs(4,1,1,2) = 1;
        coefs(4,1,2,2) = sqrt(2)/2;
        coefs(4,1,3,2) = 1;
        coefs(4,2,1,2) = 1;
        coefs(4,2,2,2) = sqrt(2)/2;
        coefs(4,2,3,2) = 1;
        
        [ GIFTmesh{indexPatch} ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW);
        
    end
    
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
    
elseif strcmp(object_type, 'cube_with_hole_C0')
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p=2;
    q=2;
    r=1;
    vol1=[];
    load('cube_with_hole_C0.mat')
    [ GIFTmesh{1} ] = genGIFTmesh3D(solidSW.knots{1}, solidSW.knots{2}, solidSW.knots{3}, solidSW.coefs, solidSW.order(1)-1, solidSW.order(2)-1, solidSW.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{2} ] = genGIFTmesh3D(solidSE.knots{1}, solidSE.knots{2}, solidSE.knots{3}, solidSE.coefs, solidSE.order(1)-1, solidSE.order(2)-1, solidSE.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{3} ] = genGIFTmesh3D(solidNW.knots{1}, solidNW.knots{2}, solidNW.knots{3}, solidNW.coefs, solidNW.order(1)-1, solidNW.order(2)-1, solidNW.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    [ GIFTmesh{4} ] = genGIFTmesh3D(solidNE.knots{1}, solidNE.knots{2}, solidNE.knots{3}, solidNE.coefs, solidNE.order(1)-1, solidNE.order(2)-1, solidNE.order(3)-1, numberElementsU, numberElementsV, numberElementsW);
    figure
    nsub=10;
    nrbplot(solidSW,[nsub,nsub,nsub])
    hold on
    nrbplot(solidSE,[nsub,nsub,nsub])
    hold on
    nrbplot(solidNW,[nsub,nsub,nsub])
    hold on
    nrbplot(solidNE,[nsub,nsub,nsub])
%    pause
    
    
elseif strcmp(object_type, 'LShaped_bracket')
    %Multipatch bracket with 18 patches
    numPatches = 18;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    numberElementsW = 1;
    p = 2;
    q = 2;
    r = 1;
    zMax = 1;
    
    a = sqrt(2)/2;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    knotW = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [0; 1; 0];
    coefs(1:3,1,3) = [-a; a; 0];
    coefs(1:3,2,1) = [2; 2; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [-2; 2; 0];
    coefs(1:3,3,1) = [4; 4; 0];
    coefs(1:3,3,2) = [0; 4; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    figure
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    
    
    GIFTmesh{1} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    
    %%%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-1; 0; 0];
    coefs(1:3,1,3) = [-a; -a; 0];
    coefs(1:3,2,1) = [-2; 2; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-2; -2; 0];
    coefs(1:3,3,1) = [-4; 4; 0];
    coefs(1:3,3,2) = [-4; 0; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    
    GIFTmesh{2} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    %%%%%%%%%%%%%%
    %third patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [0; -1; 0];
    coefs(1:3,1,3) = [a; -a; 0];
    coefs(1:3,2,1) = [-2; -2; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [2; -2; 0];
    coefs(1:3,3,1) = [-4; -4; 0];
    coefs(1:3,3,2) = [0; -4; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{3} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %Fourth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a;-a; 0];
    coefs(1:3,1,2) = [1; 0; 0];
    coefs(1:3,1,3) = [a; a; 0];
    coefs(1:3,2,1) = [2; -2; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [2; 2; 0];
    coefs(1:3,3,1) = [4; -4; 0];
    coefs(1:3,3,2) = [4; 0; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{4} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    
    %%%%%%%%%%%%%%
    %Fifth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [2; 2; 0];
    coefs(1:3,1,3) = [4; 4; 0];
    coefs(1:3,2,1) = [1; 0; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [4; 0; 0];
    coefs(1:3,3,1) = [a; -a; 0];
    coefs(1:3,3,2) = [2; -2; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8)*a; %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{5} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    
    %%%%%%%%%%%%%%
    %Sixth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; -a; 0];
    coefs(1:3,1,2) = [2; -2; 0];
    coefs(1:3,1,3) = [4; -4; 0];
    coefs(1:3,2,1) = [0; -1; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [0; -4; 0];
    coefs(1:3,3,1) = [-a; -a; 0];
    coefs(1:3,3,2) = [-2; -2; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{6} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %Seventh patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [-2; -2; 0];
    coefs(1:3,1,3) = [-4; -4; 0];
    coefs(1:3,2,1) = [-1; 0; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-4; 0; 0];
    coefs(1:3,3,1) = [-a; a; 0];
    coefs(1:3,3,2) = [-2; 2; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{7} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %Eighth patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-2; 2; 0];
    coefs(1:3,1,3) = [-4; 4; 0];
    coefs(1:3,2,1) = [0; 1; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [0; 4; 0];
    coefs(1:3,3,1) = [a; a; 0];
    coefs(1:3,3,2) = [2; 2; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{8} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %9th patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    clear coefs
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [-4;0;0];
    coefs(1:3,2,1) = [-0.853553390593274*4; 0.353553390593274*4; 0];
    coefs(1:3,3,1) = [-0.603553390593274*4; 0.603553390593274*4;0];
    coefs(1:3,1,2) = [-12;0;0];
    coefs(1:3,2,2) = [-12;8;0];
    coefs(1:3,3,2) = [-12;12;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV,knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{9} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % 10th patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    
    coefs(1:3,1,1) = [-0.603553390593274*4;0.603553390593274*4;0];
    coefs(1:3,2,1) = [-0.353553390593274*4;0.853553390593274*4;0];
    coefs(1:3,3,1) = [0;4;0];
    coefs(1:3,1,2) = [-12;12;0];
    coefs(1:3,2,2) = [-8;12;0];
    coefs(1:3,3,2) = [0;12;0];
    coefs(4,1,1) = 0.853553390593274;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV,knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{10} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %11th patch %
    %%%%%%%%%%%%%%
    clear coefs
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    a = sqrt(2)/2;
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [0; 1; 0];
    coefs(1:3,1,3) = [-a; a; 0];
    coefs(1:3,2,1) = [2; 2; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [-2; 2; 0];
    coefs(1:3,3,1) = [4; 4; 0];
    coefs(1:3,3,2) = [0; 4; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{11} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %patch 12
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-1; 0; 0];
    coefs(1:3,1,3) = [-a; -a; 0];
    coefs(1:3,2,1) = [-2; 2; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-2; -2; 0];
    coefs(1:3,3,1) = [-4; 4; 0];
    coefs(1:3,3,2) = [-4; 0; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{12} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %13th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [0; -1; 0];
    coefs(1:3,1,3) = [a; -a; 0];
    coefs(1:3,2,1) = [-2; -2; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [2; -2; 0];
    coefs(1:3,3,1) = [-4; -4; 0];
    coefs(1:3,3,2) = [0; -4; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{13} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %14th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a;-a; 0];
    coefs(1:3,1,2) = [1; 0; 0];
    coefs(1:3,1,3) = [a; a; 0];
    coefs(1:3,2,1) = [2; -2; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [2; 2; 0];
    coefs(1:3,3,1) = [4; -4; 0];
    coefs(1:3,3,2) = [4; 0; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = a;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{14} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    
    
    %%%%%%%%%%%%%%
    %15th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    
    coefs(1:3,1,1) = [a; a; 0];
    coefs(1:3,1,2) = [2; 2; 0];
    coefs(1:3,1,3) = [4; 4; 0];
    coefs(1:3,2,1) = [1; 0; 0];
    coefs(1:3,2,2) = [4; 0; 0];
    coefs(1:3,2,3) = [4; 0; 0];
    coefs(1:3,3,1) = [a; -a; 0];
    coefs(1:3,3,2) = [2; -2; 0];
    coefs(1:3,3,3) = [4; -4; 0];
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8)*a; %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{15} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %16th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:3,1,1) = [a; -a; 0];
    coefs(1:3,1,2) = [2; -2; 0];
    coefs(1:3,1,3) = [4; -4; 0];
    coefs(1:3,2,1) = [0; -1; 0];
    coefs(1:3,2,2) = [0; -4; 0];
    coefs(1:3,2,3) = [0; -4; 0];
    coefs(1:3,3,1) = [-a; -a; 0];
    coefs(1:3,3,2) = [-2; -2; 0];
    coefs(1:3,3,3) = [-4; -4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{16} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    
    %%%%%%%%%%%%%%
    %17th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; -a; 0];
    coefs(1:3,1,2) = [-2; -2; 0];
    coefs(1:3,1,3) = [-4; -4; 0];
    coefs(1:3,2,1) = [-1; 0; 0];
    coefs(1:3,2,2) = [-4; 0; 0];
    coefs(1:3,2,3) = [-4; 0; 0];
    coefs(1:3,3,1) = [-a; a; 0];
    coefs(1:3,3,2) = [-2; 2; 0];
    coefs(1:3,3,3) = [-4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (1/a-8); %correct the translation because the coefs are using projective coordinates
    %coefs(2,2,1) = a;
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{17} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
    %%%%%%%%%%%%%%
    %18th patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 2;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    
    coefs(1:3,1,1) = [-a; a; 0];
    coefs(1:3,1,2) = [-2; 2; 0];
    coefs(1:3,1,3) = [-4; 4; 0];
    coefs(1:3,2,1) = [0; 1; 0];
    coefs(1:3,2,2) = [0; 4; 0];
    coefs(1:3,2,3) = [0; 4; 0];
    coefs(1:3,3,1) = [a; a; 0];
    coefs(1:3,3,2) = [2; 2; 0];
    coefs(1:3,3,3) = [4; 4; 0];
    
    
    coefs(1,:,:,:) = coefs(1,:,:,:)-8;
    coefs(1,2,1) = (((1/a-8)*a)-1); %correct the translation because the coefs are using projective coordinates
    
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = a;
    coefs(4,2,2) = 1;
    coefs(4,2,3) = 1;
    coefs(4,3,1) = 1;
    coefs(4,3,2) = 1;
    coefs(4,3,3) = 1;
    
    coefs(:,:,:,2) = coefs(:,:,:,1);
    coefs(3,:,:,2) = zMax.*coefs(4,:,:,2);
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV knotW});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    coefs = solid2.coefs;
    
    hold on
    vol = nrbmak(coefs,{knotU, knotV, knotW});
    nrbkntplot(vol)
    GIFTmesh{18} = genGIFTmesh3D( knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW );
end




