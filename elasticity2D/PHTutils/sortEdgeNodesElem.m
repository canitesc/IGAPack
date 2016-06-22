function [ nodes, elements ] = sortEdgeNodesElem( PHTelem, edge, p, q )
%outputs the nodes and element along the given edge, sorted in the order in
%which they appear in the parameter space
%edge encoding: 1 - down, 2 - right, 3 - up, 4 - left

numEdges = length(edge);
elements = cell(1, numEdges);
nodes = cell(1, numEdges);

%initialize array to store the corners of the edge elements and the edge element indices,
%for sorting purposes
verticesDown = [];
verticesRight = [];
verticesUp = [];
verticesLeft = [];

elementsDown = [];
elementsRight = [];
elementsUp = [];
elementsLeft = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);


for elemIndex=1:length(PHTelem)
    if (isempty(PHTelem(elemIndex).children))
        for edgeIndex = 1:numEdges
            
            %down edge
            if (edge(edgeIndex) == 1) && (isempty(PHTelem(elemIndex).neighbor_down))
                %add the element to the sorted list
                elementsDown = [elementsDown, elemIndex];
                %store the left-lower corner in the verticesDown array
                %(xmin)
                verticesDown = [verticesDown, PHTelem(elemIndex).vertex(1)];                
            end
            
            %right edge
            if (edge(edgeIndex) == 2) && (isempty(PHTelem(elemIndex).neighbor_right))
                %add the element to the sorted list
                elementsRight = [elementsRight, elemIndex];
                %store the right-lower corner in the verticesRight array
                %(ymin)
                verticesRight = [verticesRight, PHTelem(elemIndex).vertex(2)];                
            end
            
            %up edge
            if (edge(edgeIndex) == 3) && (isempty(PHTelem(elemIndex).neighbor_up))
                %add the element to the sorted list
                elementsUp = [elementsUp, elemIndex];
                %store the left-upper corner in the verticesUp array
                %(xmin)
                verticesUp = [verticesUp, PHTelem(elemIndex).vertex(1)];                
            end
            
            %left edge
            if (edge(edgeIndex) == 4) && (isempty(PHTelem(elemIndex).neighbor_left))
                %add the element to the sorted list
                elementsLeft = [elementsLeft, elemIndex];
                %store the left-lower corner in the verticesLeft array
                %(ymin)
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
        elementsDown = elementsDown(indexSort);
        elements{edgeIndex} = elementsDown;
        for elemIndex = 1:length(verticesDown)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsDown(elemIndex)).nodesGlobal(down_nodes)];
        end
    end
    
    if (edge(edgeIndex) == 2)
        [verticesRight, indexSort] = sort(verticesRight);
        %rearrange the elements according to the sorting oder
        elementsRight=elementsRight(indexSort);
        elements{edgeIndex} = elementsRight;
        for elemIndex = 1:length(verticesRight)                                  
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsRight(elemIndex)).nodesGlobal(right_nodes)];
        end
    end
    
    if (edge(edgeIndex) == 3)
        [verticesUp, indexSort] = sort(verticesUp);        
        %rearrange the elements according to the sorting oder
        elementsUp = elementsUp(indexSort);
        elements{edgeIndex} = elementsUp;
        for elemIndex = 1:length(verticesUp)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsUp(elemIndex)).nodesGlobal(up_nodes)];
        end
    end
    
    if (edge(edgeIndex) == 4)
        [verticesLeft, indexSort] = sort(verticesLeft);
        %rearrange the elements according to the sorting oder
        elementsLeft = elementsLeft(indexSort);
        elements{edgeIndex} = elementsLeft;
        for elemIndex = 1:length(verticesLeft)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsLeft(elemIndex)).nodesGlobal(left_nodes)];
        end
    end
    
    %remove duplicate nodes
    nodes{edgeIndex} = unique_stab(nodes{edgeIndex});
    %nodes{edgeIndex} = unique(nodes{edgeIndex},'stable');
end
