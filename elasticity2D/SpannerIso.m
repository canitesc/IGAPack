%script Spanner.m
% implement the open-end spanner problem with PHT splines
% Use Coons parametrization

close all
clear

p = 3;
q = 3;

target_rel_error = 1e-1;
targetScale = 1e-2;

addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')


%Material properties
Emod = 1e5;
nu = 0.3;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

%the traction at infinity
tx = -10;

numPatches = 6;
[PHTelem, controlPts, dimBasis, quadList] = initPHTmeshIso('spanner',numPatches,p,q);
numberElementsU = 3;
numberElementsV = 1;

[vertices, vertex2patch, patch2vertex] = genVertex2Patch2D(PHTelem, controlPts, p, q);
[edge_list] = genEdgeList(patch2vertex);

tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
    figure
    plotPHTMeshIsoMP( PHTelem, controlPts, p, q )
    toc
    
    [ PHTelem, controlPts, dimBasis, quadList ] = checkConformingIso( PHTelem, controlPts, dimBasis, edge_list, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConformingIso( PHTelem, dimBasis, vertex2patch, edge_list, p, q);
    sizeBasis
    toc
    
        
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysIsoMP( PHTelem, controlPts, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc

    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuSpannerIso(stiff, rhs, PHTelem, controlPts, p, q, tx);
    size(stiff)
    
    toc    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    toc        
    disp('Estimating the error...')       
    [quadRef, estErrorGlobTotal] = recoverDerivEstGalIsoMPAll(PHTelem, controlPts, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal       
    toc    
    
%     %uniform refinement
%     for patchIndex=1:numPatches
%         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
%     end

     disp('Plotting the solution...')
 
     plotSolPHTElasticIsoMPVM( sol0, PHTelem, controlPts, p, q, Cmat )

     toc
%     
%     disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
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
            [quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, dimBasis(patchIndex)] = refineMeshGradedIso(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, p, q, dimBasis(patchIndex));
        end    
    end
    %stop refinment if the sum of keep_ref is zero. 
    keep_refining=sum(keep_ref);                
    %keep_refining = 1;
    toc
end