%script Blade.m
% implement a turbine blade with PHT splines
% use GIFT mapping

close all
clear all

p = 3;
q = 3;
r = 3;
tic
target_rel_error = 1e-2;
addpath ./PHTutils
addpath ./ExampleData
addpath ../nurbs/inst


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
[PHTelem{1}, dimBasis, octupleList{1}] = initPHTmesh3Dgen(p,q,r,10,20,1);


GIFTmesh{1} = init3DGeometryGIFT('blade',L,W,H);
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
    
    
%     %impose boundary conditions
%     toc
%     disp('Imposing boundary conditions...')
%     [ stiff, rhs, bcdof, bcval ] = imposeDirichletHSC1(stiff, rhs, PHTelem{1}, p, q, r, bound_disp);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    % use direct solver for this problem due to singularities in the parametrization
    sol0 = stiff\rhs;
    toc
    vtuFile = ['BladeSol',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile, 1e-3)
    toc
    %
    %     toc
    %
    %     %  figure
    disp('Estimating the error...')
    [octupleRef, estErrorGlob, estError] = recoverDerivEstGal(PHTelem{1}, GIFTmesh{1}, sol0, target_rel_error, octupleList{1}, p, q, r, Cmat);
    estErrorGlob = 1;
    %
    %     toc
        %uniform refinement
        octupleRef = 1:size(octupleList{1},1);
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


