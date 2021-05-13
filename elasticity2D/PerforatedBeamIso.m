% script PerforatedBeamIso.m
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

[PHTelem, controlPts, dimBasis, quadList] = initPHTmeshNURBSIso(nrbPatches);
[vertices, vertex2patch, patch2vertex] = genVertex2Patch2D(PHTelem, controlPts, p, q);
[edge_list] = genEdgeList(patch2vertex);

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
    
    [ PHTelem, controlPts, dimBasis, quadList ] = checkConformingIso( PHTelem,...
        controlPts, dimBasis, edge_list, p, q, quadList );
    [ PHTelem, sizeBasis ] = zipConformingIso( PHTelem, dimBasis, vertex2patch,...
        edge_list, p, q);
    
    figure
    plotPHTMeshIsoMP( PHTelem, controlPts, p, q )
    axis equal
    hold off
    drawnow
    toc
    
    %sizeBasis
    toc
    
    %assemble the linear system
    disp('Assembling the linear system...')
    [ stiff, rhs ] = assembleGalerkinSysIsoMP( PHTelem, controlPts, sizeBasis, p, q, Cmat);
    
    %impose boundary conditions
    toc
    disp('Imposing boundary conditions...')
    rhs = imposeNeumannIso(rhs, PHTelem, controlPts, neumann, traction_fun, p, q);
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
    
    [quadRef, estErrorGlobTotal] = recoverDerivEstGalIsoMPDorfler(PHTelem,...
        controlPts, sol0, quadList, p, q, Cmat, targetScale);
    disp(['Estimated error is ', num2str(estErrorGlobTotal)])
    toc
    %
    
    %     %uniform refinement
    %     for patchIndex=1:numPatches
    %         quadRef{patchIndex} = 1:size(quadList{patchIndex},1);
    %     end
    
    disp('Plotting the solution...')
    plotSolPHTElasticIsoMP( sol0, PHTelem, controlPts, p, q, Cmat )
    
    toc
    
    %adaptive refinement
    indexQuad = cell(1, numPatches);
    keep_ref = ones(1, numPatches);
    for patchIndex = 1:numPatches
        indexQuad{patchIndex} = find(quadRef{patchIndex} > 0);
        
        if isempty(indexQuad{patchIndex}) || (estErrorGlobTotal < target_rel_error)
            disp(['Done refining in geometric patch ', num2str(patchIndex),...
                ' after ',num2str(num_steps), ' steps!'])
            keep_ref(patchIndex) = 0;
        else
            numNewPatches = length(indexQuad{patchIndex});
            toc
            disp(['In geometric patch ', num2str(patchIndex), ' refining ',...
                num2str(numNewPatches), ' quadruplets out of ',...
                num2str(length(quadList{patchIndex}))])
            [quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex},...
                dimBasis(patchIndex)] = refineMeshGradedIso(quadRef{patchIndex},...
                quadList{patchIndex}, PHTelem{patchIndex}, controlPts{patchIndex},...
                p, q, dimBasis(patchIndex));
        end
    end
    %stop refinment if the sum of keep_ref is zero.
    keep_refining=sum(keep_ref);
end