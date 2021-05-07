% script ElasticRectangleMPIso.m
% implement a rectangular beam with PHT splines
% use isogeometric mapping
% variable number of patches

close all
clear

p = 3;
q = 3;

numPatches = 3;

target_rel_error = 1e-3;
targetScale = .5;
addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod = 3e4;
nu = 0;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

bound_disp = 0.1;

quadList = cell(numPatches,1);
[PHTelem, controlPts, dimBasis] = initPHTmeshIso('rectangle',numPatches,p,q);
for i=1:numPatches
    quadList{i} = 2:5;    
end

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
    plotPHTMeshIsoMP(PHTelem, controlPts, p, q)
    
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
        [quadRef{patchIndex}, estErrorGlob(patchIndex)] = recoverDerivEstGalIsoMP(PHTelem{patchIndex}, controlPts{patchIndex}, sol0, target_rel_error/numPatches, quadList{patchIndex}, p, q, Cmat, targetScale);
    end
    
    estErrorGlobTotal = sum(estErrorGlob)
    
    toc
        
%     %uniform refinement
%     for patchIndex=1:numPatches
%         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
%     end

    disp('Plotting the solution...')
    plotSolPHTElasticIsoMP( sol0, PHTelem, controlPts, p, q, Cmat )
    
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
            [quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, dimBasis(patchIndex)] = refineMeshIso(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex}, p, q, dimBasis(patchIndex));
        end    
    end
    
    %stop refinment if the sum of keep_ref is zero. 
    keep_refining=sum(keep_ref);
end


