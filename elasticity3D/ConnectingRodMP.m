%script ConnectingRod.m
% implement a connecting rod as in Phu's Nitsche coupling code with PHT splines
% use Multipatch
% use GIFT mapping

close all
clear all

p = 3;
q = 3;
r = 3;

numPatches = 25;

tic
target_rel_error = 1e-3;
targetScale = 0.5;
resultsArray = [];
addpath ./PHTutils
addpath ./ExampleData


%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

E0           = 4e5;  % Young's modulus
nu0          = 0.3;  % Poisson's ratio
bound_press = 1;   %imposed displacement on the boundary

% Define the elasticity (compliance) matrix
Cmat=zeros(6,6);
Cmat(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
Cmat(4:6,4:6)=E0/(1+nu0)*eye(3)/2;


tic
dimBasis = zeros(1, numPatches);
GIFTmesh = init3DGeometryGIFTMP('ConnectingRod',L,W,H,numPatches);

octupleList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh3D(p,q,r);
    octupleList{i} = 2:9;
end

patchBoundaries = {1, 2, 6, 5; 2, 3, 6, 5; 1, 4, 2, 4; 3, 6, 2, 4; [4,2,6],5,[6,2,5],[5,4,6]; 4,7,2,4; 6,9,2,4; [7,5,9],8,[6,2,5],[5,4,6]; 5,10,3,4;...
    10,11,2,3;11,16,5,6;11,21,6,5;16,17,2,4;21,22,2,4;[17,11,22],12,[6,2,5],[5,4,6];17,18,2,4;22,23,2,4;[18,12,23],13,[6,2,5],[5,4,6];...
    16,20,4,2;21,25,4,2;[20,11,25],15,[6,4,5],[5,2,6];[20,18],19,[4,2],[2,4];[25,23],24,[4,2],[2,4];[19,15,24,13],14,[6,4,5,2],[5,2,6,4]};


keep_refining = 1;
num_steps = 0;
%
while keep_refining
    num_steps = num_steps + 1;
    toc
    datetime
    disp(['Step ', num2str(num_steps)])
    
    % plotPHTMesh3DMP(PHTelem, GIFTmesh)
    
    [ PHTelem, dimBasis, octupleList ] = checkConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r, octupleList );
    
    [ PHTelem, sizeBasis ] = zipConforming3D( PHTelem, dimBasis, patchBoundaries, p, q, r);
    
    
    sizeBasis
    if sizeBasis>6e5
        error('Reached 1.8M DOFS... good bye!')
    end
    
    %
    toc
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, r, Cmat, octupleList);
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    [ stiff, rhs, bcdof, bcval ] = imposeDirichletConRodMP(stiff, rhs, PHTelem, GIFTmesh, p, q, r, bound_press);
    
    size(stiff)
    
    toc
    disp('Solving the linear system...')
    % sol0 = stiff\rhs;
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    toc
    stiff = M1*stiff*M1;
    
    
    %alpha = max(sum(abs(stiff),2)./diag(stiff))-2
    alpha=10
    L1 = ichol(stiff, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    toc
    [sol0,fl1,rr1,it1,rv1] = pcg(stiff,M1*rhs,1e-14,num_steps*1000,L1,L1');
    sol0 = M1*sol0;
    fprintf('PCG exited with flag %d\n', fl1)
    fprintf('Residual value: %1.15g\n', rr1)
    fprintf('Number of iterations: %d\n', it1)
    
    toc
    vtuFile = ['ConnectingRodSolutionNoGap20_p=',num2str(p),'_step',num2str(num_steps),'.vtu'];
    plotStressDisp3DVM_20pt(PHTelem, GIFTmesh, sol0, p, q, r, Cmat, vtuFile,0)
    toc
    
    [octupleRef,estErrorGlobTotal]=recoverDerivEstGalMPDorfler(PHTelem, GIFTmesh, sol0, octupleList, p, q, r, Cmat, targetScale);
    estErrorGlobTotal
    resultsArray = [resultsArray; p, length(sol0), estErrorGlobTotal]
    resultFile = ['ConnectingRodResult_p=',num2str(p),'_step',num2str(num_steps),'.mat'];
    save(resultFile, 'resultsArray')
    
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


