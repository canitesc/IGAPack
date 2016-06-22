function [patchList, PHTelem, dimBasis] = refineMesh(patchRef, patchList, PHTelem, p, q, dimBasis)
%refines the patches marked by patchRef

numElem = length(PHTelem);
indexPatch = find(patchRef > 0);
numNewPatches = length(indexPatch);
newPatchList = zeros(size(patchList,1)+numNewPatches,4);


%sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(patchList,1));
for i=1:size(patchList,1)
    levelList(i) = PHTelem(patchList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
patchList = patchList(sortedIndex, :);
patchRef = patchRef(sortedIndex);

patchCounter = 0;
for i=1:length(patchRef)
    if patchRef(i)
        curPatchElem = patchList(i,:);        
        warningFlag = 0;
        for elemIndex=patchList(i,:);
            [ ~,  ~, ~, ~, ~, ~, ~, warningFlagTemp ] = checkNeighbors( PHTelem, elemIndex );
            if warningFlagTemp
                warningFlag = 1;
            end                
        end
        if warningFlag
            disp(['Skipping refinement in quad ', num2str(i)])
            continue
        end
        
        [PHTelem, dimBasis] = crossInsert(PHTelem, curPatchElem, dimBasis, p, q );
        newPatchList(patchCounter+1,:) = numElem+1:numElem+4;
        newPatchList(patchCounter+2,:) = numElem+5:numElem+8;
        newPatchList(patchCounter+3,:) = numElem+9:numElem+12;
        newPatchList(patchCounter+4,:) = numElem+13:numElem+16;
        numElem = numElem + 16;
        patchCounter = patchCounter + 4;
    else
        newPatchList(patchCounter+1,:) = patchList(i,:);
        patchCounter = patchCounter + 1;
    end                    
end
 patchList = newPatchList;
