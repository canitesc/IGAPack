function [newBasisSet, dimBasis ] = newBasisIndices1D( elmIn, p, dimBasis )
%returns the entries of element nodes after adding new basis vertices

[ west, center, east ] = getCornerIndices1D( p );

alpha = floor((p-1)/2);

%copy the entries in elmIn to each element in newBasisSet
newBasisSet = cell(1,2);
for i=1:2
    newBasisSet{i} = elmIn;
end

basisCounter = dimBasis;

%add the deleted basis functions
newBasisSet{1}(center) = elmIn(center);
%add the other new basis functions
newBasisSet{1}(east) = basisCounter+1:basisCounter+(alpha+1);
newBasisSet{2}(west) = basisCounter+1:basisCounter+(alpha+1);
basisCounter = basisCounter + (alpha+1);
newBasisSet{2}(center) =  basisCounter+1:basisCounter+(p-2*alpha-1);
basisCounter = basisCounter + (p-2*alpha-1);

%update dimension of basis space
dimBasis = basisCounter;

