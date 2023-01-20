function [ PHTelem, dimBasis, tupleList] = initPHTmesh1D( p,numElemU)
%initialize the PHT geometry on coarse mesh, with numElemU elements

alpha = floor((p-1)/2);

knotU = [zeros(1,p+1), (1:numElemU)/numElemU, ones(1,p)];

%repeat the interior knots p-alpha times
rep_knotU = linspace(0,1,numElemU+1);
rep_knotU = rep_knotU(2:end-1);
rep_knotU = repmat(rep_knotU,1,p-alpha-1);

knotU = sort([knotU, rep_knotU]);

[C_u, ~] = bezierExtraction(knotU,p);

%compute the number of basis functions in each direction
lenU = length(knotU)-p-1;

dimBasis = lenU;
numElements = numElemU;

PHTelem = struct;
%initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];


%loop through each element and compute the element-node connectivities
elementCounter = 0;
tupleList = zeros(numElements, 2);

for i=1:length(knotU)-1
    if (knotU(i+1)>knotU(i))  %the knotspan has non-zero area
        elementCounter = elementCounter + 1;
        PHTelem(elementCounter).parent = [];
        PHTelem(elementCounter).children = [];
        PHTelem(elementCounter).vertex = [knotU(i), knotU(i+1)];
        
        %we add the nodes from i-p...i in the u direction        
        currow = (i-p):i;
                
        PHTelem(elementCounter).nodes=currow;
        PHTelem(elementCounter).level = 0;
        tupleList(elementCounter,:) = numElements + (2*(elementCounter-1)+1:2*elementCounter);
    end
end


%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = 1:numElements;

for i=1:numElemU
    elementIndex = indexMatrix(i);
    PHTelem(elementIndex).C = C_u(:,:,i);
    
    if i>1
        PHTelem(elementIndex).neighbor_left = indexMatrix(i-1);
    end
    if i<numElemU
        PHTelem(elementIndex).neighbor_right = indexMatrix(i+1);
    end
    
end
[ PHTelem, dimBasis ] = crossInsert1D( PHTelem, 1:numElements, dimBasis, p );

