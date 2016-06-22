%script CubeWithHole.m
% implement a cube with hole with PHT splines
% use C1 GIFT mapping

close all
clear all

p = 3;
q = 3;
r = 3;
tic
target_rel_error = 1e-3;
targetScale = 0.5;
addpath ./PHTutils
addpath ./ExampleData


%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

E0           = 1e3;  % Young's modulus
nu0          = 0.3;  % Poisson's ratio
bound_trac = 1;   %imposed traction on the boundary

% Define the elasticity (compliance) matrix
Cmat=zeros(6,6);
Cmat(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
Cmat(4:6,4:6)=E0/(1+nu0)*eye(3)/2;


tic
disp('Initializing the PHT geometry...')
[PHTelem{1}, dimBasis,octupleList{1}] = initPHTmesh3Dgen(p,q,r,2,2,1);


GIFTmesh{1} = init3DGeometryGIFT('cube_with_hole',L,W,H);

toc


keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    disp(['Step ', num2str(num_steps)])
    [ PHTelem, dimBasis, octupleList ] = checkConforming3D( PHTelem, dimBasis, [], p, q, r, octupleList );

    plotPHTMesh3DMP(PHTelem, GIFTmesh)
    %
    toc
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, dimBasis, p, q, r, Cmat, octupleList);
    
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeNeumannDirichletCubeWithHole(stiff, rhs, PHTelem{1}, GIFTmesh{1}, p, q, r, bound_trac);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
  %  sol0 = stiff\rhs;
    
    alpha = max(sum(abs(stiff),2)./diag(stiff))-2
    L1 = ichol(stiff, struct('type','ict','droptol',1e-3,'diagcomp',alpha));          
    [sol0,fl1,rr1,it1,rv1] = pcg(stiff,rhs,1e-14,num_steps*1000,L1,L1');            
    fprintf('PCG exited with flag %d\n', fl1)
    fprintf('Residual value: %1.15g\n', rr1)
    fprintf('Number of iterations: %d\n', it1)
    clear stiff L1
    toc
    vtuFile = ['CubeWithHoleSol',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
    toc
    
    disp('Estimating the error...')
    [octupleRef,estErrorGlob]=recoverDerivEstGalMPDorfler(PHTelem, GIFTmesh, sol0, octupleList, p, q, r, Cmat, targetScale);
    estErrorGlob
    %
        toc
        %uniform refinement
  %      octupleRef = 1:size(octupleList,1);
    %     disp('Plotting the solution...')
    %     plotSolPHTElastic( sol0, PHTelem, GIFTmesh, p, q, Cmat )
    %     %
    %     pause
    %     %calculate the error norms
        toc
        [l2relerr, h1relerr]=calcErrorNormsPHT3DCubeHole( sol0, PHTelem{1}, GIFTmesh{1}, p, q,r, Cmat );
        %l2relerr
        h1relerr
        disp(['Effectivity index: ',num2str(estErrorGlob/h1relerr)])
    %
    %adaptive refinement
    indexOctuple = find(octupleRef{1} > 0);
    
    if isempty(indexOctuple) || (estErrorGlob < target_rel_error)
        disp(['Done refining after ',num2str(num_steps), ' steps!'])
        keep_refining = 0;
    else
        numNewPatches = length(indexOctuple);
        toc
        disp(['Refining ',num2str(numNewPatches), ' octuples out of ', num2str(size(octupleList{1},1))])
        [octupleList{1}, PHTelem{1}, dimBasis] = refineMesh3D(octupleRef{1}, octupleList{1}, PHTelem{1}, p, q, r, dimBasis);
    end
    toc
end


