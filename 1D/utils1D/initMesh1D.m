function [ IGAelem, dimBasis, tupleList] = initMesh1D( p,numElemU, alpha)
%initialize the IGA mesh, with numElemU elements, alpha continuity



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

IGAelem = struct;
%initialize the neighbor connectivity lists
IGAelem.neighbor_left = [];
IGAelem.neighbor_right = [];


%loop through each element and compute the element-node connectivities
elementCounter = 0;
tupleList = zeros(numElements, 2);

for i=1:length(knotU)-1
    if (knotU(i+1)>knotU(i))  %the knotspan has non-zero area
        elementCounter = elementCounter + 1;
        IGAelem(elementCounter).parent = [];
        IGAelem(elementCounter).children = [];
        IGAelem(elementCounter).vertex = [knotU(i), knotU(i+1)];
        
        %we add the nodes from i-p...i in the u direction
        currow = (i-p):i;
        
        IGAelem(elementCounter).nodes=currow;
        IGAelem(elementCounter).level = 0;
    end
end


%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = 1:numElements;

for i=1:numElemU
    elementIndex = indexMatrix(i);
    IGAelem(elementIndex).C = C_u(:,:,i);
    
    if i>1
        IGAelem(elementIndex).neighbor_left = indexMatrix(i-1);
    end
    if i<numElemU
        IGAelem(elementIndex).neighbor_right = indexMatrix(i+1);
    end
    
end

