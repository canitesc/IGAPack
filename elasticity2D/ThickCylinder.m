%script PlateHole.m
% implement the plate with a hole problem with PHT splines
% Use C0 2-patch GIFT mapping

close all
clear

p = 3;
q = 3;

target_rel_error = 1e-6;
%targetScale = 1 --> uniform refinement
%tagetScale < 1 --> more graded refinement
targetScale = 0.75;

addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod = 1e5;
nu = 0.3;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

%boundary pressure
bound_press = 10;

%define radius of the hole and side length
rad = 1;
L = 4;

numPatches = 1;
GIFTmesh = init2DGeometryGIFTMP('thick_cylinder',1,1,1);
dimBasis = zeros(1, numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
    %[ PHTelem{i}, dimBasis(i), quadList{i}] = initPHTmeshGen( p,q, 1, 2 )
end

%patches 1 and 2 are connected at the right(2)-left(4) edges respectively
patchBoundaries = {};

tic
resultsArray = [];
keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
  %      figure
  %  plotPHTMeshMP( PHTelem, GIFTmesh )
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, patchBoundaries, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, patchBoundaries, p, q);
    sizeBasis
    toc
    figure
    plotPHTMeshMP( PHTelem, GIFTmesh )
%     plotRecovPatch
%     toc
%     pause
        
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc

    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletNeuThickCylinder(stiff, rhs, PHTelem{1}, GIFTmesh{1}, p, q, bound_press);
    size(stiff)
    
    toc    
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
    
    toc        
    disp('Estimating the error...')
    %quadList{1} = 19:22
    %[quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaAll3(PHTelem, GIFTmesh, sol0, target_rel_error, quadList, p, q, Cmat, targetScale);
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    
    
    estErrorGlobTotal
    
    toc    
%     
%     disp('Plotting the error...')
%     
%     plotErrorPHTElasticPH( sol0, PHTelem, GIFTmesh, p, q, Cmat, rad, tx )
%     %
%     toc
%     
    %calculate the error norms        
    disp('Calculating error norms...')
    rad_int = 1;
    rad_ext = 4;
    [l2relerr, h1relerr] = calcErrorNormsThickCylinder( sol0, PHTelem{1}, GIFTmesh{1}, p, q, Cmat, Emod, nu, rad_int, rad_ext, bound_press );

    l2relerr
    h1relerr

    disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
%    pause
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
    %resultsArray = [resultsArray; size(stiff,1), l2relerr, h1relerr, estErrorGlobTotal]
    toc
end