function [ elements ] = sortEdgeElem( PHTelem, edge)
%outputs the nodes and element along the given edge, sorted in the order in
%which they appear in the parameter space
%edge encoding: 1 - down, 2 - right, 3 - up, 4 - left

numEdges = length(edge);
elements = cell(1, numEdges);

%initialize array to store the corners of the edge elements and the edge element indices,
%for sorting purposes
verticesDown = [];
verticesRight = [];
verticesUp = [];
verticesLeft = [];


for elemIndex=1:length(PHTelem)
    if (isempty(PHTelem(elemIndex).children))
        for edgeIndex = 1:numEdges
            
            %down edge
            if (edge(edgeIndex) == 1) && (isempty(PHTelem(elemIndex).neighbor_down))
                %add the element to the sorted list
                elements{edgeIndex} = [elements{edgeIndex}, elemIndex];
                %store the left-lower corner in the verticesDown array
                verticesDown = [verticesDown, PHTelem(elemIndex).vertex(1)];
            end
            
            %right edge
            if (edge(edgeIndex) == 2) && (isempty(PHTelem(elemIndex).neighbor_right))
                %add the element to the sorted list
                elements{edgeIndex} = [elements{edgeIndex}, elemIndex];
                %store the right-lower corner in the verticesDown array
                verticesRight = [verticesRight, PHTelem(elemIndex).vertex(2)];                
            end
            
            %up edge
            if (edge(edgeIndex) == 3) && (isempty(PHTelem(elemIndex).neighbor_up))
                %add the element to the sorted list
                elements{edgeIndex} = [elements{edgeIndex}, elemIndex];
                %store the left-upper corner in the verticesDown array
                verticesUp = [verticesUp, PHTelem(elemIndex).vertex(1)];
            end
            
            %left edge
            if (edge(edgeIndex) == 4) && (isempty(PHTelem(elemIndex).neighbor_left))
                %add the element to the sorted list
                elements{edgeIndex} = [elements{edgeIndex}, elemIndex];
                %store the left-upper corner in the verticesDown array
                verticesLeft = [verticesLeft, PHTelem(elemIndex).vertex(2)];                
            end
        end
    end
end

%sort the elements according to vertex location
for edgeIndex=1:numEdges
    if (edge(edgeIndex) == 1)
        [verticesDown, indexSort] = sort(verticesDown);
        %rearrange the elements according to the sorting oder
        elements{edgeIndex} = elements{edgeIndex}(indexSort);        
    end
    
    if (edge(edgeIndex) == 2)
        [verticesRight, indexSort] = sort(verticesRight);
        %rearrange the elements according to the sorting oder
        elements{edgeIndex} = elements{edgeIndex}(indexSort);        
    end
    
    if (edge(edgeIndex) == 3)
        [verticesUp, indexSort] = sort(verticesUp);
        %rearrange the elements according to the sorting oder
        elements{edgeIndex} = elements{edgeIndex}(indexSort);        
    end
    
    if (edge(edgeIndex) == 4)
        [verticesLeft, indexSort] = sort(verticesLeft);
        %rearrange the elements according to the sorting oder
        elements{edgeIndex} = elements{edgeIndex}(indexSort);        
    end    
end





