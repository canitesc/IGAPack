function [ GIFTmesh] = init2DGeometryGIFTMP(object_type,L,W,numPatches)
%creates a 2d GIFT mesh and associated control points, knotvectors, Bezier
%extraction operators, etc.


if strcmp(object_type, 'rectangle')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,numPatches+1);
    
    for patchIndex = 1:numPatches
        
        %set the dimensions of the patch
        patchMinX = xVertices(patchIndex);
        patchMaxX = xVertices(patchIndex+1);
        patchMinY = 0;
        patchMaxY = W;
        
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
        coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
        coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
        coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
        coefs(4,1,1) = 1;
        coefs(4,1,2) = 1;
        coefs(4,2,1) = 1;
        coefs(4,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        
        GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
        
    end
elseif strcmp(object_type, 'plate_hole')
    
    %initialize two-patch geometry
    GIFTmesh = cell(2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = rad.*[-1;0;0];
    coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
    coefs(1:3,3,1) = rad.*[-0.603553390593274; 0.603553390593274; 0];
    coefs(1:3,1,2) = side_fac.*[-4;0;0];
    coefs(1:3,2,2) = side_fac.*[-4;2;0];
    coefs(1:3,3,2) = side_fac.*[-4;4;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 2;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = rad.*[-0.603553390593274;0.603553390593274;0];
    coefs(1:3,2,1) = rad.*[-0.353553390593274;0.853553390593274;0];
    coefs(1:3,3,1) = rad.*[0;1;0];
    coefs(1:3,1,2) = side_fac.*[-4;4;0];
    coefs(1:3,2,2) = side_fac.*[-2;4;0];
    coefs(1:3,3,2) = side_fac.*[0;4;0];
    coefs(4,1,1) = 0.853553390593274;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    GIFTmesh{2} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
elseif strcmp(object_type, 'plate_holeC1')
    numberElementsU = 2;
    numberElementsV = 1;
    p = 2;
    q = 2;
    knotU = [0, 0, 0, 0.5, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    weights = [1,(1+1/sqrt(2))/2,(1+1/sqrt(2))/2,1, 1,1,1,1, 1,1,1,1]';
    controlPoints = [-1, 0; -1, sqrt(2)-1; 1-sqrt(2), 1; 0, 1; -2.5, 0; -2.5, 0.75; -0.75, 2.5; 0, 2.5; -4, 0; -4, 4; -4, 4; 0, 4];
    %the number of control points in the u and v directions
    lenU = length(knotU)-3;  %number of basis functions in the  u direction
    lenV = length(knotV)-3;  %number of basis functions in the v direction
    coordinates = [controlPoints, weights];
    
    %convert control points/weights to coefs format in NURBS toolbox
    coefs = zeros(4,lenU,lenV);
    for i=1:lenU
        for j=1:lenV
            index = (j-1)*lenU+i;
            xcoord = coordinates(index,1);
            ycoord = coordinates(index,2);
            wght = coordinates(index,3);
            coefs(1,i,j) = xcoord*wght;
            coefs(2,i,j) = ycoord*wght;
            coefs(4,i,j) = wght;
        end
    end
    GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
elseif strcmp(object_type, 'thick_cylinder')    
    p=1;
    q=2;
    numberElementsU = 1;
    numberElementsV = 1;
    rad_int = 1;
    rad_ext = 4;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = rad_int*[1;0;0];
    coefs(1:3,1,2) = rad_int*[sqrt(2)/2; sqrt(2)/2; 0];
    coefs(1:3,1,3) = rad_int*[0;1;0];
    coefs(1:3,2,1) = rad_ext*[1;0;0];
    coefs(1:3,2,2) = rad_ext*[sqrt(2)/2;sqrt(2)/2;0];
    coefs(1:3,2,3) = rad_ext*[0;1;0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = sqrt(2)/2;
    coefs(4,1,3) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = sqrt(2)/2;
    coefs(4,2,3) = 1;
    
    knotU = [0 0 1 1];
    knotV =  [0 0 0 1 1 1];  
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'Lshape_ex')
    %L-shaped domain with exact solution discretized with 3 patches
    numPatches = 3;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [-1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,1) = [1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,2) = [0; sqrt(2); 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1/sqrt(2); 1/sqrt(2); 0];
    coefs(1:3,2,1) = [1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,2) = [sqrt(2); 0; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %third patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %define geometry on coarsest mesh
    knotU = [0, 0, 1, 1];
    knotV = [0, 0, 1, 1];
    
    coefs(1:3,1,1) = [0; 0; 0];
    coefs(1:3,1,2) = [1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,1) = [-1/sqrt(2); -1/sqrt(2); 0];
    coefs(1:3,2,2) = [0; -sqrt(2); 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'Lshape_exC1')
    numberElementsU = 2;
    numberElementsV = 1;
    p = 2;
    q = 2;
    knotU = [0, 0, 0, 0.5, 1, 1, 1];
    knotV = [0, 0, 0, 1, 1, 1];
    
    coefs(1:4,1,1) = [-1;1;0;1];
    coefs(1:4,2,1) = [-1;-1;0;1];
    coefs(1:4,3,1) = [-1;-1;0;1];
    coefs(1:4,4,1) = [1;-1;0;1];
    coefs(1:4,1,2) = [-0.65; 1; 0; 1];
    coefs(1:4,2,2) = [-0.7; 0; 0; 1];
    coefs(1:4,3,2) = [0; -0.7; 0; 1];
    coefs(1:4,4,2) = [1; -0.65; 0; 1];
    coefs(1:4,1,3) = [0; 1; 0; 1];
    coefs(1:4,2,3) = [0; 0; 0; 1];
    coefs(1:4,3,3) = [0; 0; 0; 1];
    coefs(1:4,4,3) = [1; 0; 0; 1];
    
    surf = nrbmak(coefs, {knotU, knotV});
    rotMat = vecrotz(3*pi/4);
    surf = nrbtform(surf, rotMat);
    coefs = surf.coefs;
    
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
elseif strcmp(object_type, 'CantileverBeam')
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %patch 1
    
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = W;
    patchMinY = 2*W;
    patchMaxY = 3*W;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 2
    
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = W;
    patchMinY = W;
    patchMaxY = 2*W;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 3
    
    %set the dimensions of the patch
    patchMinX = W;
    patchMaxX = W+L;
    patchMinY = W;
    patchMaxY = 2*W;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 3
    
    %set the dimensions of the patch
    patchMinX = W;
    patchMaxX = W+L;
    patchMinY = W;
    patchMaxY = 2*W;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 4
    
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = W;
    patchMinY = 0;
    patchMaxY = W;
    
    %initialize geometry on coarsest mesh
    coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
    coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
    coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
    coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
    coefs(4,1,1) = 1;
    coefs(4,1,2) = 1;
    coefs(4,2,1) = 1;
    coefs(4,2,2) = 1;
    
    knotU = [0 0 1 1];
    knotV = [0 0 1 1];
    
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'EdgeCrackMultiPatch')
    
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,3);
    yVertices = linspace(0,W,3);
    
    patchCounter = 0;
    
    for patchIndexY = 1:(numPatches-2)
        for patchIndexX = 1:(numPatches-2)
            %set the dimensions of the patch
            patchCounter = patchCounter + 1;
            patchMinX = xVertices(patchIndexX);
            patchMaxX = xVertices(patchIndexX+1);
            patchMinY = yVertices(patchIndexY);
            patchMaxY = yVertices(patchIndexY+1);
            
            %initialize geometry on coarsest mesh
            coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
            coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
            coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
            coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{patchCounter} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            
            
        end
    end
    
    tempGIFTmesh = GIFTmesh{4};
    GIFTmesh{4} = GIFTmesh{3};
    GIFTmesh{3} = tempGIFTmesh;
elseif strcmp(object_type, 'spanner')
    p=2;
    q=2;
    numberElementsU = 3;
    numberElementsV = 1;
    load('spanner.mat')
    knotU = surf.knots{1};
    knotV = [0,0,0,1,1,1];
    scale_fac = 516/12;
    surf.coefs(1:3,:,:) = scale_fac*surf.coefs(1:3,:,:);
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,1:3), p, q, numberElementsU, numberElementsV );
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,3:5), p, q, numberElementsU, numberElementsV );
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,5:7), p, q, numberElementsU, numberElementsV );
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,7:9), p, q, numberElementsU, numberElementsV );
    GIFTmesh{5} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,9:11), p, q, numberElementsU, numberElementsV );
    GIFTmesh{6} = genGIFTmesh( knotU, knotV, surf.coefs(:,:,11:13), p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'circle')
    p=2;
    q=2;
    numberElementsU = 1;
    numberElementsV = 1;
    knotU = [0,0,0,1,1,1];
    knotV = [0,0,0,1,1,1];
    coefs = [-sqrt(2)/4, -sqrt(2)/2, -sqrt(2)/4, 0, 0, 0, sqrt(2)/4, sqrt(2)/2, sqrt(2)/4; sqrt(2)/4, 0, -sqrt(2)/4, sqrt(2)/2, 0, -sqrt(2)/2, sqrt(2)/4, 0, -sqrt(2)/4; zeros(1,9);...
        1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1];
    coefs(1,:) = coefs(1,:).*coefs(4,:);
    coefs(2,:) = coefs(2,:).*coefs(4,:);
    coefs = reshape(coefs,[4,3,3]);
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type,'LShaped_bracket')
    %Multipatch bracket with 4 patches
    numPatches = 10;
    GIFTmesh = cell(numPatches,1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
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
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
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
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
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
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
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
    
    GIFTmesh{5} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    
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
    
    
    GIFTmesh{6} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    
    GIFTmesh{7} = genGIFTmesh(knotU,  knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %Eigth patch %
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
    
    GIFTmesh{8} = genGIFTmesh(knotU,knotV,coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %9th patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    coefs(1:3,3,1) = [-0.603553390593274*4;0.603553390593274*4;0];
    coefs(1:3,1,2) = [-12;0;0];
    coefs(1:3,2,2) = [-12;8;0];
    coefs(1:3,3,2) = [-12;12;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{9} =  genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    trans = vectrans([-12.0 -8.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{10} =  genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%
    %11th patch %
    %%%%%%%%%%%%%%
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
    
    %coefs(1,:,:,:) = coefs(1,:,:,:)-20;
    %coefs(2,:,:,:) = coefs(2,:,:,:)-12;
    %coefs(1,1,2) = -20; %correct the translation because the coefs are using projective coordinates
    %coefs(1,2,2) = (1/a-12);
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{11} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{12} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
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
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{13} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
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
    
    trans = vectrans([-20.0 -12.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{14} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV  );
    
    
    
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
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{15} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    
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
    
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{16} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{17} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
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
    
    trans = vectrans([-12.0 -20.0 0.0]);
    solid1 = nrbmak(coefs,{knotU knotV});
    % then, obtain the remaining by transformation
    solid2 = nrbtform(solid1, trans);
    
    GIFTmesh{18} = genGIFTmesh( knotU, knotV, solid2.coefs, p, q, numberElementsU, numberElementsV );
    
    
end
