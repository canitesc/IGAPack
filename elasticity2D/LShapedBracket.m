%script LShapedBracket.m
% implement the L-Shape domain with holes using PHT splines
% uses multipatches


close all
clear variables

p = 5;
q = 5;

target_rel_error = 1e-3;
targetScale = 1;

addpath ./PHTutils
addpath ./example_data
addpath ../nurbs/inst


%Material properties
Emod=1e5;
nu=0.3;

bound_disp = 2;

%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
%Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
a=0.5;

numPatches = 18;
GIFTmesh = init2DGeometryGIFTMP('LShaped_bracket');
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

%patchBoundaries cell array
%format for each row: {source patch(es), target patch, source patch edge(s), target patch edge(s)}
%patchBoundaries = {1, 2, 3, 1; 2, 3, 3, 1; [3,1], 4, [3, 1], [1,3]};
patchBoundaries = {1,4,1,3;
    4,3,1,3;
    [1,3],2,[3,1],[1,3];
    2, 5, 2, 3;
    5, 6, 2, 4;
    5, 8, 4, 2;
    [6,8], 7, [2,4], [4,2];
    7, 10, 3, 2;
    10, 9, 4, 2;
    9, 11, 4, 2;
    11, 14, 1, 3;
    11, 12, 3, 1;
    [12,14], 13, [3,1], [1,3];
    13, 18, 2, 3;
    18, 15, 2, 4;
    18, 17, 4, 2;
    [15,17], 16, [2,4], [4,2];
    };
tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])   
    toc        
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    
    figure
    plotPHTMeshMP( PHTelem, GIFTmesh )
    hold off
    sizeBasis
    toc
    
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs] = imposeDirichletNeuLSBracketMP(stiff, rhs, PHTelem, p, q, bound_disp);
    
    
    toc    
    disp('Solving the linear system...')
    size(stiff)
    sol0 = stiff\rhs;
    
    toc
    disp('Estimating the error...')
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);    
    estErrorGlobTotal
    %
    toc
    disp('Plotting the solution...')
    plotSolPHTElasticMPVM( sol0, PHTelem, GIFTmesh, p, q, Cmat, 5 )
    
    toc
    
    %refinement
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