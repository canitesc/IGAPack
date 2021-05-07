%script CantileverBeam.m
% implement the cantilever beam problem with PHT splines
% implement the problem from Augararde's realisitc cantilever beam paper  doi: 10.1016/j.finel.2008.01.010
% uses multipatches


close all
clear

p = 3;
q = 3;

target_rel_error = 1e-2;
targetScale = 0.5;

addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod=1000;
nu=0.25;
bound_press = 2;

shearMod=Emod/(2*(1+nu));
kappa=3-4*nu;
%elasticity matrix (plane stress)
%Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];


L = 8;
D = 2;
numPatches = 4;
GIFTmesh = init2DGeometryGIFTMP('CantileverBeam',L,D,numPatches);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

patchBoundaries = {1, 2, 1, 3; 2, 3, 2, 4; 2, 4, 1, 3};


resultsArray = [];
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
    %sizeBasis
    toc
     figure
    plotPHTMeshMP( PHTelem, GIFTmesh )
    hold off
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs] = imposeDirichletNeuCantBeam(stiff, rhs, PHTelem, GIFTmesh, p, q, L, D, bound_press);
    
    condition_number_estimate = condest(stiff)
    
    toc
    
    disp('Solving the linear system...')
    size(stiff)
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));    
    LHS = M1*stiff*M1;    
    sol0 = LHS\(M1*rhs);
    sol0 = M1*sol0;
    toc
    
    disp('Plotting the solution...')
    plotSolPHTElasticMP( sol0, PHTelem, GIFTmesh, p, q, Cmat, 1 )
    
    disp('Estimating the error...')
    %[quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal

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
            disp(['In geometric patch ', num2str(patchIndex), ' refining ',num2str(numNewPatches), ' quadruplets out of ', num2str(size(quadList{patchIndex},1))])
            [quadList{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, p, q, dimBasis(patchIndex));
        end
    end
    
    %stop refinment if the sum of keep_ref is zero.
    keep_refining=sum(keep_ref);
    %resultsArray = [resultsArray; size(stiff,1), l2relerr, h1relerr, estErrorGlobTotal]
end