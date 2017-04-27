% script ElasticRectangleNeuMP.m
% implement a rectangular beam with PHT splines
% use GIFT mapping
% Neumann boundary conditions on the right edge

close all
clear all

p = 4;
q = 4;

numPatches = 1;

target_rel_error = 1e-5;
targetScale = 0.5;
addpath ./PHTutils
addpath ./example_data
addpath ../nurbs/inst

%Material properties
Emod = 3e4;
nu = 0.3;
%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

bound_trac = [100, 0];

L=1;
W=1;

dimBasis = zeros(1, numPatches);


GIFTmesh = init2DGeometryGIFTMP('rectangle',L,W,numPatches);

quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;    
end

patchBoundaries = cell(numPatches-1, 4);
for indexPatch = 1:numPatches-1
    %define the boundary between patch i and patch i+1
    patchBoundaries{indexPatch,1} = indexPatch;
    patchBoundaries{indexPatch,2} = indexPatch+1;
    patchBoundaries{indexPatch,3} =  2;
    patchBoundaries{indexPatch,4} = 4;
end
tic

keep_refining = 1;
num_steps = 0;

while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
    figure
    plotPHTMeshMP(PHTelem, GIFTmesh)
   % pause
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
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletESNeuMP(stiff, rhs, PHTelem, GIFTmesh, p, q, bound_trac);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    sol0 = stiff\rhs;
       
    toc
  
    %  figure
    disp('Estimating the error...')    
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    estErrorGlobTotal
    
    toc
        
%     %uniform refinement
%     for patchIndex=1:numPatches
%         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
%     end
    
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


