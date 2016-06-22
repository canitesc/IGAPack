%script LShapePanelIso.m
% implement the L-Shaped panel problem with PHT splines
% implement exact solution from R. Simpson's paper/code
% uses multipatches


close all
%clear all

p = 3;
q = 3;

target_rel_error = 1e-4;
targetScale = 0.5;

addpath ./PHTutils
addpath ./example_data
addpath ../nurbs/inst

%Material properties
Emod=1e5;
nu=0.3;

shearMod=Emod/(2*(1+nu));
kappa=3-4*nu;
%elasticity matrix (plane stress)
%Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];


numPatches = 3;
[PHTelem, controlPts, dimBasis, quadList] = initPHTmeshIso('Lshape_ex',numPatches,p,q);

%patches 1 and 2 are connected at the down(1)-left(4) edges respectively
%patches 2 and 3 are connected at the down(1)-left(4) edges respectively
patchBoundaries = {1, 2, 1, 4; 2, 3, 1, 4};



tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
%     for i=1:numPatches
%         [ PHTelem{i}] = scaleBasis(PHTelem{i}, dimBasis(i));
%     end
    disp(['Step ', num2str(num_steps)])
    figure
    plotPHTMeshIsoMP( PHTelem, controlPts, p, q )    
    hold off
    toc
    
    [ PHTelem, controlPts, dimBasis, quadList ] = checkConformingIso( PHTelem, controlPts, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    %sizeBasis
    toc
    
    %assemble the linear system
   disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysIsoMP( PHTelem, controlPts, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs] = imposeDirichletNeuLSExSymIso(stiff, rhs, PHTelem, controlPts, p, q);
    
    %condition_number_estimate = condest(stiff)
    
    toc
    
    disp('Solving the linear system...')
    condest(stiff)
    toc
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    toc
    LHS = M1*stiff*M1;
    
    toc
    condest(LHS)
    toc
    sol0 = LHS\(M1*rhs);
    toc
    sol0 = M1*sol0;
    %pause
    size(stiff)
    toc   
    disp('Estimating the error...')
    
    quadRef = cell(1, numPatches);
    estErrorGlob = zeros(1, numPatches);
    
    for patchIndex=1:numPatches
        [quadRef{patchIndex}, estErrorGlob(patchIndex)] = recoverDerivEstGalIsoMP(PHTelem{patchIndex}, controlPts{patchIndex}, sol0, target_rel_error/numPatches, quadList{patchIndex}, p, q, Cmat, targetScale);
    end
    
    estErrorGlobTotal = sum(estErrorGlob)
    %
    toc
    %
    
%     %uniform refinement
%     for patchIndex=1:numPatches
%         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
%     end

    disp('Plotting the solution...')
    plotSolPHTElasticIsoMP( sol0, PHTelem, controlPts, p, q, Cmat )
    
    toc
    %calculate the error norms
    disp('Calculating error norms...')
    [l2relerr, h1relerr] = calcErrorNormsLSIso( sol0, PHTelem, controlPts, p, q, Cmat, shearMod, kappa );
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
end