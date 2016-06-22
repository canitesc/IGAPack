%script ElasticCUBEPHT.m
% implement a cube domain with PHT splines
% use GIFT mapping

close all
clear all

p = 5;
q = 5;
r = 5;
tic
targetScale = 0.5;
target_rel_error = 1e-2;
addpath ./PHTutils
addpath ./ExampleData


%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

E0           = 1e5;  % Young's modulus
nu0          = 0.3;  % Poisson's ratio
bound_disp = 0.1;   %imposed displacement on the boundary

% Define the elasticity (compliance) matrix
Cmat=zeros(6,6);
Cmat(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
Cmat(4:6,4:6)=E0/(1+nu0)*eye(3)/2;


tic
[PHTelem{1}, dimBasis] = initPHTmesh3D(p,q,r);
toc
GIFTmesh{1} = init3DGeometryGIFT('cube',L,W,H);
octupleList{1} = 2:9;

keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    toc
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
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletEC(stiff, rhs, PHTelem{1}, p, q, r, bound_disp);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    toc
    alpha = 10;
    stiff = M1*stiff*M1';
    toc
    L1 = ichol(stiff, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    toc
    [sol0,fl1,rr1,it1,rv1] = pcg(stiff,M1*rhs,1e-14,num_steps*1000,L1,L1');
    sol0 = M1*sol0;
    fprintf('PCG exited with flag %d\n', fl1)
    fprintf('Residual value: %1.15g\n', rr1)
    fprintf('Number of iterations: %d\n', it1)
    
    toc
    vtuFile = ['ElasticCubeC1Sol',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
    toc

    % error estimation
    [octupleRef, estErrorGlob, estError] = recoverDerivEstGalMPAll(PHTelem, GIFTmesh, sol0, target_rel_error, octupleList, p, q, r, Cmat, targetScale);
    estErrorGlob
    
    %adaptive refinement
    indexOctuple = find(octupleRef{1} > 0);
    
    if isempty(indexOctuple) || (estErrorGlob < target_rel_error)
        disp(['Done refining after ',num2str(num_steps), ' steps!'])
        keep_refining = 0;
    else
        numNewPatches = length(indexOctuple);
        toc
        disp(['Refining ',num2str(numNewPatches), ' octuples out of ', num2str(size(octupleList,1))])
        [octupleList{1}, PHTelem{1}, dimBasis] = refineMesh3D(octupleRef{1}, octupleList{1}, PHTelem{1}, p, q, r, dimBasis);
    end
    toc
end


