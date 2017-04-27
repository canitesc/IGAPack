function [ PHTelem, controlPts, dimBasis, quadList ] = checkConformingIso( PHTelem, controlPts, dimBasis, patchBoundaries, p, q, quadList )
%checks that the patches are conforming and if needed makes them conforming through mesh refinement
%patchBoundaries format:
% patchA, patchB, edgeA, edgeB
% first patchA should be 1
%edge format: 1-down, 2-right, 3-up, 4-left

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHTelem);

%create/set nodesGlobal entries in all patches to be equal to local nodes
%entries
keepChecking = 1;
while keepChecking
    keepChecking = 0;
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end
    
    for boundaryIndex = 1:numBoundaries
        %get the nodes on the boundary edge in patchA and patchB
        patchA = patchBoundaries{boundaryIndex,1};
        patchB = patchBoundaries{boundaryIndex,2};
        edgeA = patchBoundaries{boundaryIndex,3};
        edgeB = patchBoundaries{boundaryIndex,4};
        quadListA = quadList{patchA};
        quadListB = quadList{patchB};
        
        [elementsA] = sortEdgeElem( PHTelem{patchA}, edgeA);
        [elementsB] = sortEdgeElem( PHTelem{patchB}, edgeB);
        
        [quadRefA, quadRefB] = makeConforming(PHTelem{patchA}, PHTelem{patchB}, elementsA{1}, elementsB{1}, edgeA, edgeB, quadListA, quadListB);
        indexQuadA = find(quadRefA > 0);
        indexQuadB = find(quadRefB > 0);
        
        if ~isempty(indexQuadA)
            keepChecking = 1;
            numNewPatches = length(indexQuadA);
            disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchB)])
            [quadList{patchA}, PHTelem{patchA}, controlPts{patchA}, dimBasis(patchA)] = refineMeshIso(quadRefA, quadListA, PHTelem{patchA}, controlPts{patchA}, p, q, dimBasis(patchA));
        end
        
        if ~isempty(indexQuadB)
            keepChecking = 1;
            numNewPatches = length(indexQuadB);
            disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewPatches), ' quadruplets to keep conformity with patch ', num2str(patchA)])
            [quadList{patchB}, PHTelem{patchB}, controlPts{patchB}, dimBasis(patchB)] = refineMeshIso(quadRefB, quadListB, PHTelem{patchB}, controlPts{patchB}, p, q, dimBasis(patchB));
        end
    end
end



