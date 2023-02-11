function [ IGAelem, controlPts, dimBasis, nrb, numberElementsUV] = initNURBSmeshIso(object_type,numPatches,p,q )
%initialize the NURBS geometry on coarse mesh

%start with num_elements_u, num_elements_v uniformly distributed knot spans
%in each parametric direction

if strcmp(object_type, 'plate_hole')
    %initialize two-patch geometry
    controlPts = cell(numPatches,1);
    IGAelem = cell(numPatches, 1);
    nrbs = cell(numPatches, 1);
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
    
    [controlPts{1}, IGAelem{1}, dimBasis(1), nrb{1}] = genControlPtsNURBS( knotU, ...
        knotV, coefs, p, q, numberElementsU, numberElementsV );
    numberElementsUV{1} = [numberElementsU, numberElementsV];

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
    
    [controlPts{2}, IGAelem{2}, dimBasis(2), nrb{2}] = genControlPtsNURBS( knotU, ...
        knotV, coefs, p, q, numberElementsU, numberElementsV );
    numberElementsUV{2} = [numberElementsU, numberElementsV];
    
end
