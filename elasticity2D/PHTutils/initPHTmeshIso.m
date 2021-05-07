function [ PHTelem, controlPts, dimBasis, quadList ] = initPHTmeshIso(object_type,numPatches,p,q )
%initialize the PHT geometry on coarse mesh

%start with num_elements_u, num_elements_v uniformly distributed knot spans
%in each parametric direction

if strcmp(object_type, 'rectangle')
    
    L = 8;
    W = 2;
    
    controlPts = cell(numPatches,1);
    PHTelem = cell(numPatches, 1);
    numberElementsU = 1;
    numberElementsV = 1;
    dimBasis = zeros(1, numPatches);
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,numPatches+1);
    quadList = cell(1,numPatches);
    for patchIndex = 1:numPatches
        
        dimBasis(patchIndex) = (p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1)
        
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
        
        [controlPts{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex), quadList{patchIndex}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
        
    end
elseif strcmp(object_type, 'plate_hole')
    
    %initialize two-patch geometry
    controlPts = cell(numPatches,1);
    PHTelem = cell(numPatches, 1);
    dimBasis = ones(1, numPatches)*(p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %data for the first patch%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    
    %initialize geometry on coarsest mesh
    %define geometry on coarsest mesh
    knotU = [0, 0, 0, 1, 1, 1];
    knotV = [0, 0, 1, 1];
    
    %initialize geometry on coarsest mesh
    rad = 1;
    side_fac = 1;
    coefs(1:3,1,1) = rad.*[-1;0;0];
    coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
    coefs(1:3,3,1) = rad.*[-0.603553390593274;0.603553390593274;0];
    coefs(1:3,1,2) = side_fac.*[-4;0;0];
    coefs(1:3,2,2) = side_fac.*[-4;2;0];
    coefs(1:3,3,2) = side_fac.*[-4;4;0];
    coefs(4,1,1) = 1;
    coefs(4,2,1) = 0.853553390593274;
    coefs(4,3,1) = 0.853553390593274;
    coefs(4,1,2) = 1;
    coefs(4,2,2) = 1;
    coefs(4,3,2) = 1;
    
    %GIFTmesh{1} =  genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    [controlPts{1}, PHTelem{1}, dimBasis(1), quadList{1}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for second patch%
    %%%%%%%%%%%%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    
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
    
    [controlPts{2}, PHTelem{2}, dimBasis(2), quadList{2}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'Lshape_ex')
    %initialize three-patch geometry
    numPatches = 3;
    controlPts = cell(numPatches,1);
    PHTelem = cell(numPatches, 1);
    dimBasis = ones(1, numPatches)*(p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1);
    
    %%%%%%%%%%%%%
    %first patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    
    %initialize geometry on coarsest mesh
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
    
    [controlPts{1}, PHTelem{1}, dimBasis(1), quadList{1}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %do the initial refinement   
    
    
    %%%%%%%%%%%%%
    %second patch%
    %%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    
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
    
    [controlPts{2}, PHTelem{2}, dimBasis(2), quadList{2}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    %do the initial refinement
    
    %%%%%%%%%%%%%%
    %third patch %
    %%%%%%%%%%%%%%
    numberElementsU = 1;
    numberElementsV = 1;
    
    
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
    
    [controlPts{3}, PHTelem{3}, dimBasis(3), quadList{3}] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
elseif strcmp(object_type, 'spanner')
    load('spanner.mat')
    knotU = surf.knots{1};
    knotV = [0,0,0,1,1,1];
    scale_fac = 516/12;
    numberElementsU = 3;
    numberElementsV = 1;
    surf.coefs(1:3,:,:) = scale_fac*surf.coefs(1:3,:,:);
    [controlPts{1}, PHTelem{1}, dimBasis(1), quadList{1}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,1:3), p, q, numberElementsU, numberElementsV );
    [controlPts{2}, PHTelem{2}, dimBasis(2), quadList{2}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,3:5), p, q, numberElementsU, numberElementsV );
    [controlPts{3}, PHTelem{3}, dimBasis(3), quadList{3}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,5:7), p, q, numberElementsU, numberElementsV );
    [controlPts{4}, PHTelem{4}, dimBasis(4), quadList{4}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,7:9), p, q, numberElementsU, numberElementsV );
    [controlPts{5}, PHTelem{5}, dimBasis(5), quadList{5}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,9:11), p, q, numberElementsU, numberElementsV );
    [controlPts{6}, PHTelem{6}, dimBasis(6), quadList{6}] = genControlPtsNoRep( knotU, knotV, surf.coefs(:,:,11:13), p, q, numberElementsU, numberElementsV );    
   
end
% dimBasis

% dimBasis
% pause
