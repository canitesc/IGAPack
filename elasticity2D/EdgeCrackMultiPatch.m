% script EdgeCrackMultiPatch.m
% implement a Edge crack with PHT splines
% use GIFT mapping
% use 4 patches

close all
clear

p = 3;
q = 3;

numPatches = 4;

target_rel_error = 1e-3;
targetScale = 0.5;

addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod = 1e3;
nu = 0.3;
%elasticity matrix (plane stress)
%Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
%elasticity matrix (plane strain)
stressState = 'PLANE_STRAIN';
Cmat = Emod/((1+nu)*(1-2*nu))*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];


force = 1;

L=1;
W=2;

dimBasis = zeros(1, numPatches);


GIFTmesh = init2DGeometryGIFTMP('EdgeCrackMultiPatch',L,W,numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;    
end

patchBoundaries = {1, 2, 2, 4; 2, 3, 3, 1; 3, 4, 4, 2};

tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])

    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);        
    sizeBasis
    figure
    plotPHTMeshMP(PHTelem, GIFTmesh)
    
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletECMP(stiff, rhs, PHTelem, GIFTmesh, p, q, force);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));    
    LHS = M1*stiff*M1;    
    sol0 = LHS\(M1*rhs);
    sol0 = M1*sol0;
       
    toc
  
    %  figure
    disp('Estimating the error...')    
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal
    
    toc
        
    disp('Plotting the solution...')
    plotSolPHTElasticMPVM( sol0, PHTelem, GIFTmesh, p, q, Cmat, 100 )
    
   %vtuFile = 'EdgeCrackMultiPatch.vtu';
   %PlotStressDispMP 
     
   %calculate the error norms
    disp('Calculating error norms...')
    [l2relerr, h1relerr] = calcErrorNormsECMP( sol0, PHTelem, GIFTmesh, p, q, Cmat, force, Emod, nu, stressState);
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
end


