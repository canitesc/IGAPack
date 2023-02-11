%script PlateHoleIso.m
% implement the plate with a hole problem with NURBS
% Use C0 2-patch isoparametric mapping

close all
clear

p = 3;
q = 3;

target_rel_error = 1e-5;
resultsArray =[];


addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod = 1e5;
nu = 0.3;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

%the traction at infinity
tx = 10;

%define radius of the hole and side length
rad = 1;
L = 4;

numPatches = 2;

[IGAelem, controlPts, dimBasis, nurbs, numberElementsUV] = initNURBSmeshIso('plate_hole', numPatches, p, q);
[vertices, vertex2patch, patch2vertex] = genVertex2Patch2D(IGAelem, controlPts, p, q);
[edge_list] = genEdgeList(patch2vertex);


tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc    
    disp(['Step ', num2str(num_steps)])
    close all
    figure
   
    toc  
    [ IGAelem, sizeBasis ] = zipConforming( IGAelem, dimBasis, vertex2patch, ...
        edge_list, p, q);
    sizeBasis
    toc
     plotPHTMeshIsoMP( IGAelem, controlPts, p, q )
     %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysIsoMP( IGAelem, controlPts, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPIso(stiff, rhs, IGAelem, ...
        controlPts, p, q, rad, tx);
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    toc
    disp('Plotting the solution...')
    plotSolPHTElasticIsoMP( sol0, IGAelem, controlPts, p, q, Cmat )
    
    %calculate the error norms        
    disp('Calculating error norms...')
    [l2relerr, h1relerr] = calcErrorNormsPHIso( sol0, IGAelem, controlPts, ...
        p, q, Cmat, Emod, nu, rad, tx );
    l2relerr
    h1relerr
    
    % uniform refinement    
    for patchIndex = 1:numPatches
        [nurbs{patchIndex}] = refineMeshIsoUniformNURBS(nurbs{patchIndex});
        knotU = nurbs{patchIndex}.knots{1};
        knotV = nurbs{patchIndex}.knots{2};
        coefs = nurbs{patchIndex}.coefs;
        numberElementsU = 2*numberElementsUV{patchIndex}(1);
        numberElementsV = 2*numberElementsUV{patchIndex}(2);
        numberElementsUV{patchIndex} = [numberElementsU, numberElementsV];
        [ controlPts{patchIndex}, IGAelem{patchIndex}, dimBasis(patchIndex), nurbs{patchIndex}] = ...
            genControlPtsNURBS( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );

    end
    
    % stop refining if the error is less than the target
    if h1relerr<target_rel_error
        keep_refining = 0;
    end
    toc
end