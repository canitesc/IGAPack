function [ nodes, elements ] = sortEdgeNodesElem( PHTelem, edge, p, q , flagDir)
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
        [verticesDown, indexSort] = sort(flagDir(edgeIndex)*verticesDown);
        %rearrange the elements according to the sorting oder
        elementsDown = elementsDown(indexSort);
        elements{edgeIndex} = elementsDown;
        if flagDir(edgeIndex)==-1
            down_nodes_dir = flip(down_nodes);
        else
            down_nodes_dir = down_nodes;
        end
        for elemIndex = 1:length(verticesDown)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsDown(elemIndex)).nodes(down_nodes_dir)];
        end
    end
    
    if (edge(edgeIndex) == 2)
        [verticesRight, indexSort] = sort(flagDir(edgeIndex)*verticesRight);
        %rearrange the elements according to the sorting oder
        elementsRight=elementsRight(indexSort);
        elements{edgeIndex} = elementsRight;
        if flagDir(edgeIndex)==-1
            right_nodes_dir = flip(right_nodes);
        else
            right_nodes_dir = right_nodes;
        end
        for elemIndex = 1:length(verticesRight)                                  
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsRight(elemIndex)).nodes(right_nodes_dir)];
        end
    end
    
    if (edge(edgeIndex) == 3)
        [verticesUp, indexSort] = sort(flagDir(edgeIndex)*verticesUp);        
        %rearrange the elements according to the sorting oder
        elementsUp = elementsUp(indexSort);
        elements{edgeIndex} = elementsUp;
        if flagDir(edgeIndex)==-1
            up_nodes_dir = flip(up_nodes);
        else
            up_nodes_dir = up_nodes;
        end
        for elemIndex = 1:length(verticesUp)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsUp(elemIndex)).nodes(up_nodes_dir)];
        end
    end
    
    if (edge(edgeIndex) == 4)
        [verticesLeft, indexSort] = sort(flagDir(edgeIndex)*verticesLeft);
        %rearrange the elements according to the sorting oder
        elementsLeft = elementsLeft(indexSort);
        elements{edgeIndex} = elementsLeft;
        if flagDir(edgeIndex)==-1
            left_nodes_dir = flip(left_nodes);
        else
            left_nodes_dir = left_nodes;
        end
        for elemIndex = 1:length(verticesLeft)
            nodes{edgeIndex} = [nodes{edgeIndex}, PHTelem(elementsLeft(elemIndex)).nodes(left_nodes_dir)];
        end
    end
    
    %remove duplicate nodes
    nodes{edgeIndex} = unique_stab(nodes{edgeIndex});
    %nodes{edgeIndex} = unique(nodes{edgeIndex},'stable');
end
