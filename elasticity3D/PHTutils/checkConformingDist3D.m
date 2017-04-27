function checkConformingDist3D( PHTelem, GIFTmesh, patchBoundaries, p)
%checks that the patches are conforming and if needed makes them conforming through mesh refinement
%additionally checks that the parametric order matches by evaluating the
%distance between the joined faces at the Gauss points
%patchBoundaries format:
% patchA, patchB, faceA, faceB
% patchA should be 1
%edge format: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up

numBoundaries = size(patchBoundaries,1);
numPatches = length(PHTelem);

keepChecking = 1;
while keepChecking
    keepChecking = 0;
    %create/set nodesGlobal entries in all patches to be equal to local nodes
    %entries
    for patchIndex = 1:numPatches
        for elemIndex = 1:length(PHTelem{patchIndex})
            PHTelem{patchIndex}(elemIndex).nodesGlobal = PHTelem{patchIndex}(elemIndex).nodes;
        end
    end
    
    for boundaryIndex = 1:numBoundaries
        %get the nodes on the boundary edge in patchA and patchB
        patchAList = patchBoundaries{boundaryIndex,1};
        patchB = patchBoundaries{boundaryIndex,2};
        faceAList = patchBoundaries{boundaryIndex,3};
        faceBList = patchBoundaries{boundaryIndex,4};
        
        for indexPatch=1:length(patchAList)
            patchA = patchAList(indexPatch);
            faceA = faceAList(indexPatch);
            faceB = faceBList(indexPatch);            
            
            [elementsA] = sortFaceElem( PHTelem{patchA}, faceA);
            [elementsB] = sortFaceElem( PHTelem{patchB}, faceB);
            boundaryIndex
            [distance] = makeConformingDist3D(PHTelem{patchA}, PHTelem{patchB}, GIFTmesh{patchA}, GIFTmesh{patchB}, elementsA, elementsB, faceA, faceB, p)

        end
    end
end



