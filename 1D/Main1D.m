%script Main1D.m
% solve -u''(x) + a1(x)*u'(x) + a0(x)*u(x) = f(x), x \in (0, L), with Dirichlet boundary conditions
% u(0) = a, u(1) = b;
%
%pre-processing

clear all
close all

addpath ./utils1D
addpath ./ExampleData
addpath ../nurbs/inst

GeoMesh = init1DGeometryBeam;


numberElements = 2; %initial number of non-empty knot-spans; this should be an even number!

target_rel_error = 1e-10;
targetScale = 0.5;

p = 4; %polynomial order


numPatches = length(GeoMesh);
PHTelem = cell(1,numPatches);
tupleList = cell(1,numPatches);
[ PHTelem{1}, dimBasis, tupleList{1}] = initPHTmesh1D( p,numberElements);
numberElements = 2*numberElements; %after cross refinement we have twice the number of elements
tupleRef = cell(1, numPatches);


keep_refining = 1;

stepCounter = 0;
tic
while (keep_refining == 1)
    stepCounter = stepCounter + 1;
    
    disp(['Step ', num2str(stepCounter)])
    
    %plotBasisPHTMesh1D(PHTelem, GeoMesh, p)
    
    %pause
    
    %calculate stiffness matrix and force vector using Bezier extraction
    numberElementsOld = numberElements;
    [lhs, rhs, numberElements] = assembleGalerkin1D(PHTelem, GeoMesh, dimBasis, p);
    
    %impose boundary conditions
    [lhs, rhs] = boundaryCond1D(lhs, rhs, PHTelem, GeoMesh, 'DN');
    
    disp('Solving linear system...')
    sol0 = lhs\rhs;
    
    disp('Estimating error...')
    
    estErrorGlob = zeros(1, numPatches);
    
    for patchIndex = 1:numPatches
        %[tupleRef{patchIndex}, estErrorGlob(patchIndex)] = recoveryGreedy1D(PHTelem{patchIndex}, GeoMesh{patchIndex}, sol0, target_rel_error, tupleList{patchIndex}, p );
        [tupleRef{patchIndex}, estErrorGlob(patchIndex)] = recoveryDoerfler1D(PHTelem{patchIndex}, GeoMesh{patchIndex}, sol0, tupleList{patchIndex}, p, targetScale);
    end
    
    estErrorGlobTotal = sum(estErrorGlob);
    
    if estErrorGlobTotal < target_rel_error
        keep_refining = 0;
        disp('Estimated error is below target threshold')
    end
    % estErr
    disp(['Estimated H1 error: ', num2str(estErrorGlobTotal)])
    if numberElements > numberElementsOld
        l2normold = l2relerr;
        h1normold = h1relerr;
        
        [l2relerr,h1relerr]=calcError1D(PHTelem,GeoMesh,sol0,p);
        
        disp('Error in square norms')
        disp(sprintf('%1.15f  %1.15f', [l2relerr,h1relerr]))
        
        oconvl2 = log(l2normold/l2relerr)/log(dimBasis/dimBasisOld);
        oconvh1 = log(h1normold/h1relerr)/log(dimBasis/dimBasisOld);
        
        disp(['Number of Elements: ', num2str(numberElements)])
        disp(['O(conv) L2: ', num2str(oconvl2), ' O(conv) H1: ', num2str(oconvh1)])
    else
        [l2relerr,h1relerr]=calcError1D(PHTelem,GeoMesh,sol0,p);
        
        disp('Error in square norms')
        disp(sprintf('%1.15f  %1.15f', [l2relerr,h1relerr]))
    end
    
    disp(['Effectivity index: ',num2str(estErrorGlobTotal/h1relerr)])
    errorArray(stepCounter,:) = [dimBasis, h1relerr, estErrorGlobTotal];
    dimBasisOld = dimBasis;
    
    disp('Plotting the actual error...')
    figure
    plotError1D( PHTelem, GeoMesh, sol0, p )

    %adaptive refinement
    indexTuple = find(tupleRef{1} > 0);
    
    if isempty(indexTuple) || (~keep_refining)
        keep_refining = 0;
        disp('Finished refining!')
    else
        numNewTuples = length(indexTuple);
        numPatches = size(tupleList,1);
        disp(['Refining ',num2str(numNewTuples), ' element tuples out of ', num2str(size(tupleList{patchIndex},1))])
        [tupleList{1}, PHTelem{1}, dimBasis] = refineMesh1D( PHTelem{1}, tupleList{1}, tupleRef{1}, p, dimBasis );
    end
    numberElementsOld = numberElements;
    toc
    
end
%pause

figure
loglog(errorArray(:,1), errorArray(:,2), '-r', errorArray(:,1), errorArray(:,3), '-b')
legend('Exact error', 'Estimated error')
title('Convergence of H1 seminorm error vs. number of degrees of freedom')
disp(['Finished after ',num2str(stepCounter), ' steps'])
%plotErrors1d( numuniqelem,elementNodes,b_net,knotVector,p,displacements )
%plotErrors1d_col( numuniqelem,elementNodes,b_net,knotVector,p,displacements, sample_points, element_int);


