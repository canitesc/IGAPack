%script Spanner.m
% implement the open-end spanner problem with PHT splines
% Use Coons parametrization

close all
clear

p = 4;
q = 4;

target_rel_error = 1e-1;
targetScale = 0.9;

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
GIFTmesh = init2DGeometryGIFTMP('spanner');
dimBasis = zeros(1, numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
numberElementsU = 3;
numberElementsV = 1;
for i=1:numPatches
    [PHTelem{i}, dimBasis(i), quadList{i}] = initPHTmeshGen(p,q,numberElementsU,numberElementsV);
end

%patches 1 and 2 are connected at the right(2)-left(4) edges respectively
patchBoundaries = {1,2,3,1;2,3,3,1;3,4,3,1;4,5,3,1;5,6,3,1};

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
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    sizeBasis
    toc
    
    
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc
    
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuSpanner(stiff, rhs, PHTelem, GIFTmesh, p, q, tx);
    size(stiff)
    toc
    
    disp('Solving the linear system...')
    size(stiff)
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    LHS = M1*stiff*M1;
    sol0 = LHS\(M1*rhs);
    sol0 = M1*sol0;
    toc
    
    disp('Estimating the error...')
    %[quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal
    toc
    
    %     %uniform refinement
    %     for patchIndex=1:numPatches
    %         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
    %     end
    
    disp('Plotting the solution...')
    plotSolPHTElasticMPVM( sol0, PHTelem, GIFTmesh, p, q, Cmat )
    
    toc
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