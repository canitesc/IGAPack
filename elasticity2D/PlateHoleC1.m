%script PlateHole.m
% implement the plate with a hole problem with PHT splines
% Use C1 1-patch GIFT mapping

close all
clear

p = 3;
q = 3;

target_rel_error = 1e-5;
targetScale = 1;

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

numPatches = 1;
GIFTmesh = init2DGeometryGIFTMP('plate_holeC1');
dimBasis = zeros(1, numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

[vertices, vertex2patch, patch2vertex] = genVertex2PatchGift2D(GIFTmesh);
[edge_list] = genEdgeList(patch2vertex);

tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
    figure
    plotPHTMeshMP( PHTelem, GIFTmesh )
    toc
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, edge_list, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, vertex2patch, edge_list, p, q);
    sizeBasis
    toc
    
        
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc

    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuPH(stiff, rhs, PHTelem{1}, GIFTmesh{1}, p, q, rad, tx);
    size(stiff)
    
    toc    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    toc        
    disp('Estimating the error...')
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    %[quadRef,estErrorGlobTotal]=recoverDerivEstGalExactErrMPPH(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale, rad, tx);

    estErrorGlobTotal
    
    toc    
    
%     %uniform refinement
%     for patchIndex=1:numPatches
%         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
%     end

    disp('Plotting the error...')
    plotErrorPHTElasticPH( sol0, PHTelem, GIFTmesh, p, q, Cmat, rad, tx )
    %
    toc
    
    %calculate the error norms        
    disp('Calculating error norms...')
    [l2relerr, h1relerr] = calcErrorNormsPH( sol0, PHTelem, GIFTmesh, p, q, Cmat, Emod, nu, rad, tx );
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
            [quadList{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, p, q, dimBasis(patchIndex));
        end    
    end
    
    %stop refinment if the sum of keep_ref is zero. 
    keep_refining=sum(keep_ref);                
    %keep_refining = 1;
    toc
end