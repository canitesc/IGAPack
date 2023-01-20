function [ tupleList, PHTelem, dimBasis] = refineMesh1D( PHTelem, tupleList, tupleRef, p, dimBasis )
%refine the PHT mesh by cross insertion in the marked elements

numElem = length(PHTelem);
indexTuple = find(tupleRef > 0);
numNewTuples = length(indexTuple);
newTupleList = zeros(size(tupleList,1)+numNewTuples,2);


%sort tupleList by level of elements to prevent coarse->fine level errors
levelList = zeros(1,size(tupleList,1));
for i=1:size(tupleList,1)
    levelList(i) = PHTelem(tupleList(i,1)).level;
end
    
[~,sortedIndex] = sort(levelList);
tupleList = tupleList(sortedIndex, :);
tupleRef = tupleRef(sortedIndex);

tupleCounter = 0;
for i=1:length(tupleRef)
    if tupleRef(i)
        curTupleElem = tupleList(i,:);        
        
        [PHTelem, dimBasis] = crossInsert1D(PHTelem, curTupleElem, dimBasis, p);
        newTupleList(tupleCounter+1,:) = numElem+1:numElem+2;
        newTupleList(tupleCounter+2,:) = numElem+3:numElem+4;      
        numElem = numElem + 4;
        tupleCounter = tupleCounter + 2;
    else
        newTupleList(tupleCounter+1,:) = tupleList(i,:);
        tupleCounter = tupleCounter + 1;
    end                    
end
tupleList = newTupleList;
