function [patchList,PHTelem,controlPts,dimBasis,markRef] = refineMeshGradedIso(patchRef,patchList,PHTelem,controlPts,p,q,dimBasis)
%refines the patches marked by patchRef
%also refine the neighbor patches if the level_neighbor<level_refined_quad

numElem = length(PHTelem);
markRef = zeros(length(PHTelem),1);
indexPatch = find(patchRef > 0);
numNewPatches = length(indexPatch);
newPatchList = zeros(size(patchList,1)+numNewPatches,4);

keepChecking = 1;
% Check the refinement level for the neighbor quad and mark it for refinement if needed
extraPatchCounter = 0;
while keepChecking
    keepChecking = 0;
    patchRefOrig = patchRef;
    for i=1:length(patchRefOrig)
        if patchRefOrig(i)
            curPatchElemSW = patchList(i,1);
            curPatchElemNE = patchList(i,4);
            curLevel = PHTelem(curPatchElemSW).level;
            if ~isempty(PHTelem(curPatchElemSW).neighbor_down)
                neighborSouth = PHTelem(curPatchElemSW).neighbor_down(1);
                if curLevel>PHTelem(neighborSouth).level
                    % disp('South')
                    [patchIndex, ~] = find(patchList==neighborSouth);
                    if patchRef(patchIndex)==0
                        patchRef(patchIndex) = 1;
                        extraPatchCounter = extraPatchCounter + 1;
                        keepChecking = 1;
                    end
                    %pause
                end
            end
            if ~isempty(PHTelem(curPatchElemSW).neighbor_left)
                %disp('West')
                neighborWest = PHTelem(curPatchElemSW).neighbor_left(1);
                if curLevel>PHTelem(neighborWest).level
                    [patchIndex, ~] = find(patchList==neighborWest);
                    if patchRef(patchIndex)==0
                        patchRef(patchIndex) = 1;
                        extraPatchCounter = extraPatchCounter + 1;
                        keepChecking = 1;
                    end
                    %pause
                end
            end
            if ~isempty(PHTelem(curPatchElemNE).neighbor_up)
                %disp('North')
                neighborNorth = PHTelem(curPatchElemNE).neighbor_up(1);
                if curLevel>PHTelem(neighborNorth).level
                    [patchIndex, ~] = find(patchList==neighborNorth);
                    if patchRef(patchIndex)==0
                        patchRef(patchIndex) = 1;
                        extraPatchCounter = extraPatchCounter + 1;
                        keepChecking = 1;
                    end
                    %pause
                end
            end
            if ~isempty(PHTelem(curPatchElemNE).neighbor_right)
                %disp('East')
                neighborEast = PHTelem(curPatchElemNE).neighbor_right(1);
                if curLevel>PHTelem(neighborEast).level
                    [patchIndex, ~] = find(patchList==neighborEast);
                    if patchRef(patchIndex)==0
                        patchRef(patchIndex) = 1;
                        extraPatchCounter = extraPatchCounter + 1;
                        keepChecking = 1;
                    end
                end
            end
        end
    end
    disp(['Extra patches refined: ',num2str(extraPatchCounter)]);
end

% Sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(patchList,1));
for i=1:size(patchList,1)
    levelList(i) = PHTelem(patchList(i,1)).level;
end

[~,sortedIndex] = sort(levelList);
patchList = patchList(sortedIndex,:);
patchRef = patchRef(sortedIndex);
tempPHTelem = PHTelem;
patchCounter = 0;
for i=1:length(patchRef)
    if patchRef(i)
        curPatchElem = patchList(i,:);
        markRef(curPatchElem) = 1;
        for iCount = 1:4
            curElem = curPatchElem(iCount);
            north = tempPHTelem(curElem).neighbor_up;
            south = tempPHTelem(curElem).neighbor_down;
            east = tempPHTelem(curElem).neighbor_right;
            west = tempPHTelem(curElem).neighbor_left;
           
            neighbor = [east,west,north,south];
            markRef(neighbor) = 1;
        end
        [PHTelem,controlPts,dimBasis] = crossInsertIso(PHTelem,controlPts,curPatchElem,dimBasis,p,q);
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
