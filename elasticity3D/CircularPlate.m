%script CircularPlate.m
% implement a clamped circular plate with PHT splines
% use GIFT mapping


close all
clear all

p = 5;
q = 5;
r = 5;

numPatches = 1;
targetScale = 0.75;
tic
target_rel_error = 1e-2;
addpath ./PHTutils
addpath ./ExampleData


%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

E0           = 30e6;  % Young's modulus
nu0          = 0.2;  % Poisson's ratio
bound_disp = 0.1;   %imposed displacement on the boundary

% Define the elasticity (compliance) matrix
Cmat=zeros(6,6);
Cmat(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
Cmat(4:6,4:6)=E0/(1+nu0)*eye(3)/2;

modelRho =  2.320;
numModes = 7;

tic
dimBasis = zeros(1, numPatches);
GIFTmesh = init3DGeometryGIFTMP('circular_plate',L,W,H,numPatches);

octupleList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh3D(p,q,r);
    octupleList{i} = 2:9;
end

patchBoundaries = cell(numPatches-1, 4);


keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
    
    plotPHTMesh3DMP(PHTelem, GIFTmesh)
    [ PHTelem, dimBasis, octupleList ] = checkConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r, octupleList );
    
    [ PHTelem, sizeBasis ] = zipConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r);
    
    
    sizeBasis
    
    %
    toc
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, mass, rhs ] = assembleGalerkinMassSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, r, Cmat, modelRho, octupleList );
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, mass, bcdof ] = imposeDirichletCircPlate(stiff, mass, PHTelem, p, q, r);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    activeDOF = setdiff(1:3*sizeBasis,bcdof);
    
    %     alpha = max(sum(abs(stiff),2)./diag(stiff))-2
    %     L1 = ichol(stiff, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    %
    %     [sol0,fl1,rr1,it1,rv1] = pcg(stiff,rhs,1e-14,num_steps*1000,L1,L1');
    %
    %     fprintf('PCG exited with flag %d\n', fl1)
    %     fprintf('Residual value: %1.15g\n', rr1)
    %     fprintf('Number of iterations: %d\n', it1)
    
    disp('Computing eigenvalues...')
    [V, omega] = eigs(stiff, mass, numModes, 'SM');
    %[omega] = eigs(stiff, mass,7,'SM');
    disp('omega');
    omega = diag(omega);
    omega = sort(sqrt(omega));
    fprintf('%12.8f\n',omega(1:numModes));
 
    toc
    for indexMode = numModes:-1:1
        sol0 = zeros(3*sizeBasis,1);
        sol0(activeDOF) = V(:,indexMode);
        vtuFile = ['CircularPlate_p=',num2str(p),'_mode=',num2str(indexMode),'_step=',num2str(num_steps),'.vtu'];
        plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile)
    end
    
    
    [octupleRef,estErrorGlobTotal]=recoverDerivEstGalMPDorfler(PHTelem, GIFTmesh, sol0, octupleList, p, q, r, Cmat, targetScale);
    estErrorGlobTotal
    
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


