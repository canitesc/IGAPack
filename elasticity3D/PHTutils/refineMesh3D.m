function [octupleList, PHTelem, dimBasis] = refineMesh3D(octupleRef, octupleList, PHTelem, p, q, r, dimBasis)
%refines the patches marked by patchRef

numElem = length(PHTelem);
indexOctuple = find(octupleRef > 0);
numNewOctuples = length(indexOctuple);
newOctupleList = zeros(size(octupleList,1)+numNewOctuples,8);


%sort patchList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(octupleList,1));
for i=1:size(octupleList,1)
    levelList(i) = PHTelem(octupleList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
octupleList = octupleList(sortedIndex, :);
octupleRef = octupleRef(sortedIndex);
octupleCounter = 0;
octupleIndexCounter = 0;
for i=1:length(octupleRef)
    if octupleRef(i)
        octupleIndexCounter = octupleIndexCounter + 1;
        curPatchElem = octupleList(i,:);        
        disp(['Refining octuple ',num2str(octupleIndexCounter), ' out of ', num2str(numNewOctuples), '...'])
        [PHTelem, dimBasis] = crossInsert3D(PHTelem, curPatchElem, dimBasis, p, q, r );
        newOctupleList(octupleCounter+1,:) = numElem+1:numElem+8;
        newOctupleList(octupleCounter+2,:) = numElem+9:numElem+16;
        newOctupleList(octupleCounter+3,:) = numElem+17:numElem+24;
        newOctupleList(octupleCounter+4,:) = numElem+25:numElem+32;
        newOctupleList(octupleCounter+5,:) = numElem+33:numElem+40;
        newOctupleList(octupleCounter+6,:) = numElem+41:numElem+48;
        newOctupleList(octupleCounter+7,:) = numElem+49:numElem+56;
        newOctupleList(octupleCounter+8,:) = numElem+57:numElem+64;
        numElem = numElem + 64;
        octupleCounter = octupleCounter + 8;
    else
        newOctupleList(octupleCounter+1,:) = octupleList(i,:);
        octupleCounter = octupleCounter + 1;
    end                    
end
 octupleList = newOctupleList;
