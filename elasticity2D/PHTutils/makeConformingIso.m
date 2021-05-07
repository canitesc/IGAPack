function [ quadRefA, quadRefB ] = makeConforming( PHTelemA, PHTelemB, elementsA, elementsB, edgeA, edgeB, quadListA, quadListB, flagDir )
%checks that patchA and patchB are conforming (the elements match up) and
%refines if necessary
%Note: assumes only 1 level of refinement is needed

%initialize quadRefA and quadRefB with zeros (no refinements)
quadRefA = zeros(1, size(quadListA,1));
quadRefB = zeros(1, size(quadListB,1));

%make list of coordinates of element vertices along the patch boundaries
numElementsA = length(elementsA);
numElementsB = length(elementsB);

vertexA = zeros(1, length(elementsA));
vertexB = zeros(1, length(elementsB));

if (edgeA == 1) || (edgeA == 3) %horizontal edge
    vertIndexA = 3;
else
    vertIndexA = 4;
end

if (edgeB == 1) || (edgeB == 3) %horizontal edge
    if flagDir==1
        vertIndexB = 3;
    else
        vertIndexB = 1;
    end
else
    if flagDir==1
        vertIndexB = 4;
    else
        vertIndexB = 2;
    end
end

%loop over the elements in A
for elementIndex = 1:numElementsA    
    vertexA(elementIndex) = PHTelemA(elementsA(elementIndex)).vertex(vertIndexA);
end

%loop over the elements in B
for elementIndex = 1:numElementsB
    vertexB(elementIndex) = PHTelemB(elementsB(elementIndex)).vertex(vertIndexB);
end

if flagDir==1
    newVertexB = setdiff(vertexA, vertexB);
    newVertexA = setdiff(vertexB, vertexA);
else
    newVertexB = setdiff(1-vertexA, vertexB);
    newVertexA = setdiff(1-vertexB, vertexA);
end

vertexA = [0, vertexA];
if flagDir==1
    vertexB = [0, vertexB];
else
    vertexB = flip([1, vertexB]);
end

%find the elements and quads in A that need to be refined
for vertIndex = 1:length(newVertexA)
    for elementIndex = 1:numElementsA
        if (vertexA(elementIndex+1) > newVertexA(vertIndex)) && (vertexA(elementIndex)<newVertexA(vertIndex))
            [quadIndex, ~] = find(quadListA==elementsA(elementIndex));
            quadRefA(quadIndex) = 1; %mark the corresponding quad for refinements
            break
        end
    end
end

%find the elements and quads in B that need to be refined
for vertIndex = 1:length(newVertexB)
    if flagDir==-1
        elementsB = flip(elementsB);
    end
    for elementIndex = 1:numElementsB
        if (vertexB(elementIndex+1) > newVertexB(vertIndex)) && (vertexB(elementIndex)<newVertexB(vertIndex))
            [quadIndex, ~] = find(quadListB==elementsB(elementIndex));
            quadRefB(quadIndex) = 1; %mark the corresponding quad for refinements
            break
        end
    end
end           
