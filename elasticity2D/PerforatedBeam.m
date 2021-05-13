% script PerforatedBeam.m
%   Elasticity equation where the domain is a perforated rectangle with
% corners at (0,-1) and (10,1) with 8 holes of radius 0.2 centered at
% (2.5, -0.5), (2.5, 0.5), (4, 0), (5, -0.5), (5, 0.5),
%    (6, 0), (7.5, -0.5), (7.5, 0.5)
%    The left end has zero displacements in x and y direction
%    The right hand has a parabolic traction \tau_xy = P/2(D^2/4-y^2)

close all
clear

p = 3;
q = 3;

target_rel_error = 1e-2;
targetScale = 0.5;

addpath('./PHTutils')
addpath('./example_data')
addpath('../nurbs/inst')

%Material properties
Emod=1e5;
nu=0.3;

%elasticity matrix (plane stress)
Cmat=Emod/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

nrbPatches = genPerforatedBeamGeom();
numPatches = length(nrbPatches);

%convert to GIFT mesh
GIFTmesh = cell(1, numPatches);
for i=1:numPatches
    knotU = nrbPatches{i}.knots{1};
    knotV = nrbPatches{i}.knots{2};
    coefs = nrbPatches{i}.coefs;
    pgeom = nrbPatches{i}.order(1)-1;
    qgeom = nrbPatches{i}.order(2)-1;
    numElemU = length(unique(knotU))-1;
    numElemV = length(unique(knotV))-1;
    GIFTmesh{i} = genGIFTmesh(knotU, knotV, coefs, pgeom, qgeom, numElemU, numElemV);
end

dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
for i=1:numPatches
    [PHTelem{i}, dimBasis(i)] = initPHTmesh(p,q);
    quadList{i} = 2:5;
end

[vertices, vertex2patch, patch2vertex] = genVertex2PatchGift2D(GIFTmesh);
[edge_list] = genEdgeList(patch2vertex);
tic

% boundary condition data
dirichlet = [1, 4, 1, 1; 2, 4, 1, 1];
neumann = [3, 2; 4, 2];
param.P = 10;
param.D = 2;
traction_fun = @(coords, normal) parabolicTraction(coords, normal, param);

% adaptive refinement loop
tic
keep_refining = 1;
num_steps = 0;

while keep_refining
    close all
    num_steps = num_steps + 1;
    toc
    disp(['Step ', num2str(num_steps)])
    
    toc
    [ PHTelem, dimBasis, quadList ] = checkConforming( PHTelem, dimBasis, edge_list, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, vertex2patch, edge_list, p, q);
    
    figure
    plotPHTMeshMP( PHTelem, GIFTmesh )
    axis equal
    drawnow
    hold off
    sizeBasis
    toc
    
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysGIFTMP( PHTelem, GIFTmesh, sizeBasis, p, q, Cmat);
    
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    rhs = imposeNeumann(rhs, PHTelem, GIFTmesh, neumann, traction_fun, p, q);
    [ stiff, rhs] = imposeHomogeneousDirichlet(stiff, rhs, PHTelem, dirichlet, p, q);
    
    %condition_number_estimate = condest(stiff)
    
    toc
    
    disp('Solving the linear system...')
    condest(stiff)
    toc
    M1=spdiags(1./sqrt(diag(stiff)),0,size(stiff,1),size(stiff,2));
    toc
    LHS = M1*stiff*M1;
    
    toc
    condest(LHS)
    toc
    sol0 = LHS\(M1*rhs);
    toc
    sol0 = M1*sol0;
    %pause
    size(stiff)
    toc
    disp('Estimating the error...')
    [quadRef,estErrorGlobTotal]=recoverDerivEstGalMP2alphaDorfler(PHTelem, GIFTmesh, sol0, quadList, p, q, Cmat, targetScale);
    disp(['Estimated error is ', num2str(estErrorGlobTotal)])
    %
    
    %     %uniform refinement
    %     for patchIndex=1:numPatches
    %         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
    %     end
    
    disp('Plotting the solution...')
    plotSolPHTElasticMP( sol0, PHTelem, GIFTmesh, p, q, Cmat, 10 )
    
    toc
    
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
            disp(['In geometric patch ', num2str(patchIndex), ' refining ',num2str(numNewPatches), ' quadruplets out of ', num2str(size(quadList{patchIndex},1))])
            [quadList{patchIndex}, PHTelem{patchIndex}, dimBasis(patchIndex)] = refineMesh(quadRef{patchIndex}, quadList{patchIndex}, PHTelem{patchIndex}, p, q, dimBasis(patchIndex));
        end
    end
    
    %stop refinment if the sum of keep_ref is zero.
    keep_refining=sum(keep_ref);
end