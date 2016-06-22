function [ octupleRefA, octupleRefB ] = makeConforming3D( PHTelemA, PHTelemB, elementsA, elementsB, faceA, faceB, octupleListA, octupleListB)
%checks that patchA and patchB are conforming (the elements match up) and
%refines if necessary
%Note: assumes only 1 level of refinement is needed

%initialize octupleRefA and octupleRefB with zeros (no refinements)
octupleRefA = zeros(1, size(octupleListA,1));
octupleRefB = zeros(1, size(octupleListB,1));

%make list of coordinates of element vertices along the patch boundaries
numElementsA = length(elementsA);
numElementsB = length(elementsB);

vertexA = zeros(length(elementsA), 2);
vertexB = zeros(length(elementsB), 2);
vertexMaxA = zeros(length(elementsA), 2);
vertexMaxB = zeros(length(elementsB), 2);

if (faceA == 1) || (faceA == 3) %front/back face
    vertIndexA = [1,3];
elseif (faceA == 2) || (faceA == 4) %left/right face
    vertIndexA = [2,3];
else %up/down face
    vertIndexA = [1,2];
end

if (faceB == 1) || (faceB == 3) %front/back face
    vertIndexB = [1,3];
elseif (faceB == 2) || (faceB == 4) %left/right face
    vertIndexB = [2,3];
else %up/down face
    vertIndexB = [1,2];
end

%loop over the elements in A
for elementIndex = 1:numElementsA    
    vertexA(elementIndex,:) = [PHTelemA(elementsA(elementIndex)).vertex(vertIndexA(1)), PHTelemA(elementsA(elementIndex)).vertex(vertIndexA(2))];
    vertexMaxA(elementIndex,:) = [PHTelemA(elementsA(elementIndex)).vertex(vertIndexA(1)+3), PHTelemA(elementsA(elementIndex)).vertex(vertIndexA(2)+3)];
end

%loop over the elements in B
for elementIndex = 1:numElementsB
    vertexB(elementIndex,:) = [PHTelemB(elementsB(elementIndex)).vertex(vertIndexB(1)), PHTelemB(elementsB(elementIndex)).vertex(vertIndexB(2))];
    vertexMaxB(elementIndex,:) = [PHTelemB(elementsB(elementIndex)).vertex(vertIndexB(1)+3), PHTelemB(elementsB(elementIndex)).vertex(vertIndexB(2)+3)];
end

newVertexB = setdiff(vertexA, vertexB,'rows');
newVertexA = setdiff(vertexB, vertexA,'rows');

vertexA = [0, 0; vertexA];
vertexB = [0, 0; vertexB];



toleq =0;

%find the elements and octuples in A that need to be refined
for vertIndex = 1:size(newVertexA,1)
    for elementIndex = 1:numElementsA
        if (vertexMaxA(elementIndex,1)+toleq > newVertexA(vertIndex,1)) && (vertexA(elementIndex,1)<newVertexA(vertIndex,1)+toleq) ...
                && (vertexMaxA(elementIndex,2)+toleq > newVertexA(vertIndex,2)) && (vertexA(elementIndex,2)<newVertexA(vertIndex,2)+toleq)
            [octupleIndex, ~] = find(octupleListA==elementsA(elementIndex));
            octupleRefA(octupleIndex) = 1; %mark the corresponding octuple for refinements            
            break
        end
    end
end

%find the elements and octuples in B that need to be refined
for vertIndex = 1:size(newVertexB,1)
    for elementIndex = 1:numElementsB
        if (vertexMaxB(elementIndex,1)+toleq > newVertexB(vertIndex,1)) && (vertexB(elementIndex,1)<newVertexB(vertIndex,1)+toleq) ...
                && (vertexMaxB(elementIndex,2)+toleq > newVertexB(vertIndex,2)) && (vertexB(elementIndex,2)<newVertexB(vertIndex,2)+toleq)
            [octupleIndex, ~] = find(octupleListB==elementsB(elementIndex));
            octupleRefB(octupleIndex) = 1; %mark the corresponding octuple for refinements
            break
        end
    end
end

