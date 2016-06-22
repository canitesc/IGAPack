%script Hemisphere.m
% implement an 8th a sphere domain (pressurized sphere problem) with PHT splines
% use GIFT mapping

close all
clear all

p = 4;
q = 4;
r = 4;
tic
target_rel_error = 1e-4;
addpath ./PHTutils
addpath ./ExampleData

targetScale = 0.75;
numPatches = 1;
%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

E0           = 1e3;  % Young's modulus
nu0          = 0.3;  % Poisson's ratio
bound_press = 1;   %imposed displacement on the boundary

% Define the elasticity (compliance) matrix
Cmat=zeros(6,6);
Cmat(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
Cmat(4:6,4:6)=E0/(1+nu0)*eye(3)/2;


tic
disp('Initializing the PHT geometry...')
[PHTelem{1}, dimBasis,octupleList{1}] = initPHTmesh3Dgen(p,q,r,1,1,1);
%octupleList = 2:9;

GIFTmesh{1} = init3DGeometryGIFT('hemisphere',L,W,H);

toc
%octupleList = 2:9;
%octupleList = [3:10;11:18];

patchBoundaries = [];
keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    disp(['Step ', num2str(num_steps)])
    
    plotPHTMesh3D(PHTelem{1}, GIFTmesh{1})
    [ PHTelem, dimBasis, octupleList ] = checkConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r, octupleList );
    
    [ PHTelem, sizeBasis ] = zipConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r);
    sizeBasis = dimBasis;
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, r, Cmat, octupleList);
    
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletHemi(stiff, rhs, PHTelem{1}, GIFTmesh{1}, p, q, r, bound_press);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    %  sol0 = stiff\rhs;
    
    alpha = max(sum(abs(stiff),2)./diag(stiff))-2
    L1 = ichol(stiff, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    %L1 = ichol(stiff, struct('diagcomp',alpha));
    [sol0,fl1,rr1,it1,rv1] = pcg(stiff,rhs,1e-14,num_steps*1000,L1,L1');
    fprintf('PCG exited with flag %d\n', fl1)
    fprintf('Residual value: %1.15g\n', rr1)
    fprintf('Number of iterations: %d\n', it1)
    clear stiff L1
    toc
    vtuFile = ['HemisphereSolution',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
    toc
    %
    %     toc
    %
    %     %  figure
    [octupleRef,estErrorGlobTotal]=recoverDerivEstGalMPDorfler(PHTelem, GIFTmesh, sol0, octupleList, p, q, r, Cmat, targetScale);
    estErrorGlobTotal
        
    disp('Calculating the actual error...')
    [l2relerr, h1relerr]=calcErrorNormsPHT3DHemisphere( sol0, PHTelem{1}, GIFTmesh{1}, p, q, r, Cmat);
    l2relerr
    h1relerr
    disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
    
    toc
   
   %adaptive refinement
    indexOctuple = cell(1, numPatches);
    keep_ref = ones(1, numPatches);
    for patchIndex = 1:numPatches
        indexOctuple{patchIndex} = find(octupleRef{patchIndex} > 0);
        
        if isempty(indexOctuple{patchIndex}) || (estErrorGlobTotal < target_rel_error)
            disp(['Done refining in geometric patch ', num2str(patchIndex), ' after ',num2str(num_steps), ' steps!'])
            keep_ref(patchIndex) = 0;
        else
            numNewOctuples = length(indexOctuple{patchIndex});
            toc
            disp(['In geometric patch ', num2str(patchIndex), ' refining ',num2str(numNewOctuples), ' octuples out of ', num2str(size(octupleList{patchIndex},1))])
            [octupleList{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh3D(octupleRef{patchIndex}, octupleList{patchIndex}, PHTelem{patchIndex}, p, q, r, dimBasis(patchIndex));
        end
    end
    %stop refinment if the sum of keep_ref is zero.
    keep_refining=sum(keep_ref);
    toc
end


