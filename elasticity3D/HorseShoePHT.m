%script HorseShoePHT.m
% implement a horseshoe domain with PHT splines
% use GIFT mapping

close all
clear all

p = 4;
q = 4;
r = 4;
tic
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
disp('Initializing the PHT geometry...')
[PHTelem{1}, dimBasis, octupleList{1}] = initPHTmesh3Dgen(p,q,r,5,2,1);


GIFTmesh{1} = init3DGeometryGIFT('horseshoe',L,W,H);
toc
%octupleList = 2:9;
%octupleList = [3:10;11:18];


keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    disp(['Step ', num2str(num_steps)])
    [ PHTelem, dimBasis, octupleList ] = checkConforming3D( PHTelem, dimBasis, [], p, q, r, octupleList );
    plotPHTMesh3D(PHTelem{1}, GIFTmesh{1})
    %
    toc
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFT( PHTelem{1}, GIFTmesh{1}, dimBasis, p, q, r, Cmat, octupleList{1});
    
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletHSC1(stiff, rhs, PHTelem{1}, p, q, r, bound_disp);
    
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
    clear stiff L1
    toc
    vtuFile = ['HorseShoeSol',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
    toc
    %
    %     toc
    %
    %     %  figure
    disp('Estimating the error...')
    [octupleRef, estErrorGlob, estError] = recoverDerivEstGal(PHTelem{1}, GIFTmesh{1}, sol0, target_rel_error, octupleList{1}, p, q, r, Cmat);
    estErrorGlob
    %
    %     toc
    %     %%uniform refinement
    %     %patchRef = 1:size(patchList,1);
    %     disp('Plotting the solution...')
    %     plotSolPHTElastic( sol0, PHTelem, GIFTmesh, p, q, Cmat )
    %     %
    %     pause
    %     %calculate the error norms
    %     toc
    %
    %adaptive refinement
    indexOctuple = find(octupleRef > 0);
    
    if isempty(indexOctuple) || (estErrorGlob < target_rel_error)
        disp(['Done refining after ',num2str(num_steps), ' steps!'])
        keep_refining = 0;
    else
        numNewPatches = length(indexOctuple);
        toc
        disp(['Refining ',num2str(numNewPatches), ' octuples out of ', num2str(size(octupleList{1},1))])
        [octupleList{1}, PHTelem{1}, dimBasis] = refineMesh3D(octupleRef, octupleList{1}, PHTelem{1}, p, q, r, dimBasis);
    end
    toc
end


