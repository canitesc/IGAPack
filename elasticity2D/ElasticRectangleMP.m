% script ElasticRectangleMP.m
% implement a rectangular beam with PHT splines
% use GIFT mapping
% use 2 patches
% This is a patch test of a plate in tension where the solution has linear
% displacement and constant stresses/strains

close all
clear

p = 4;
q = 4;

numPatches = 3;

target_rel_error = 1e-3;
addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')
%Material properties
Emod = 3e4;
nu = 0;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

bound_disp = 0.1;

L=1;
W=1;

dimBasis = zeros(1, numPatches);


GIFTmesh = init2DGeometryGIFTMP('rectangle',L,W,numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [ PHTelem{i}, dimBasis(i), quadList{i}] = initPHTmeshGen( p,q, 2, 2 );
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

    
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, edge_list, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, vertex2patch, edge_list, p, q);    
        figure
    plotPHTMeshMP(PHTelem, GIFTmesh)
    sizeBasis
    toc
%    pause
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
%        pause
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletESMP(stiff, rhs, PHTelem, p, q, bound_disp);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
       
    toc
  
    %  figure
    disp('Estimating the error...')
    
    quadRef = cell(1, numPatches);
    estErrorGlob = zeros(1, numPatches);
    
    for patchIndex=1:numPatches
        [quadRef{patchIndex}, estErrorGlob(patchIndex)] = recoverDerivEstGalMP(PHTelem{patchIndex}, GIFTmesh{patchIndex}, sol0, target_rel_error, quadList{patchIndex}, p, q, Cmat);
    end
    
    estErrorGlobTotal = sum(estErrorGlob)
    
    toc
        
    %uniform refinement
    for patchIndex=1:numPatches
        quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
    end

    disp('Plotting the solution...')
    plotSolPHTElasticMP( sol0, PHTelem, GIFTmesh, p, q, Cmat )
    %
    
    %calculate the error norms
    
    
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
end


