%script PlateHoleIso.m
% implement the plate with a hole problem with PHT splines
% Use C0 2-patch isoparametric mapping

close all
clear all

p = 3;
q = 3;

target_rel_error = 1e-6;
targetScale = 0.75;
resultsArray =[];


addpath ./PHTutils
addpath ./example_data
addpath ../nurbs/inst

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

[PHTelem, controlPts, dimBasis, quadList] = initPHTmeshIso('plate_hole',numPatches,p,q);

%patches 1 and 2 are connected at the right(2)-left(4) edges respectively
patchBoundaries = {1, 2, 2, 4};

tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc    
    disp(['Step ', num2str(num_steps)])
    figure
   
    toc
    [ PHTelem, controlPts, dimBasis, quadList ] = checkConformingIso( PHTelem, controlPts, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    sizeBasis
    toc
     plotPHTMeshIsoMP( PHTelem, controlPts, p, q )
     %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysIsoMP( PHTelem, controlPts, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPHMPIso(stiff, rhs, PHTelem, controlPts, p, q, rad, tx);
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
       
    toc
  
    %  figure
    disp('Estimating the error...')       
    %[quadRef, estErrorGlobTotal] = recoverDerivEstGalIsoMPAll(PHTelem, controlPts, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    [quadRef, estErrorGlobTotal] = recoverDerivEstGalIsoMPDorfler(PHTelem, controlPts, sol0, quadList, p, q, Cmat, targetScale);

    estErrorGlobTotal
    
    toc

    %disp('Plotting the solution...')
    %plotSolPHTElasticIsoMP( sol0, PHTelem, controlPts, p, q, Cmat )
    
    %calculate the error norms        
    disp('Calculating error norms...')
    [l2relerr, h1relerr] = calcErrorNormsPHIso( sol0, PHTelem, controlPts, p, q, Cmat, Emod, nu, rad, tx );
    l2relerr
    h1relerr
    disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
    
    %adaptive refinement
    indexQuad = cell(1, numPatches);
    keep_ref = ones(1, numPatches);
    for patchIndex = 1:numPatches
        indexQuad{patchIndex} = find(quadRef{patchIndex} > 0);

        if isempty(indexQuad{patchIndex}) || (estErrorGlobTotal < target_rel_error)
            disp(['Done refining in geometric patch ', num2str(patchIndex), ' after ',num2str(num_steps), ' steps!'])
            keep_ref(patchIndex) = 0;
        else
            numNewPatches = length(indexQuad{patchIndex});
            toc       
            disp(['In geometric patch ', num2str(patchIndex), ' refining ',num2str(numNewPatches), ' quadruplets out of ', num2str(length(quadList{patchIndex}))])
            [quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, dimBasis(patchIndex)] = refineMeshIso(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, p, q, dimBasis(patchIndex));
        end    
    end
    
    %stop refinment if the sum of keep_ref is zero. 
    keep_refining=sum(keep_ref);
    %resultsArray = [resultsArray; size(stiff,1), l2relerr, h1relerr, estErrorGlobTotal]

    toc
end