function [ nodes, elements ] = sortFaceNodesElem( PHTelem, face, p, q, r )
%outputs the nodes and element along the given edge, sorted in the order in
%which they appear in the parameter space
%Note: here face is scalar (single value)
%face encoding: 1 - front, 2 - right, 3 - back, 4 - left, 5 - down, 6 - up

%initialize array to store the corners of the edge elements and the edge element indices,
%for sorting purposes

verticesFront = [];
verticesRight = [];
verticesBack = [];
verticesLeft = [];
verticesDown = [];
verticesUp = [];

elementsFront = [];
elementsRight = [];
elementsBack = [];
elementsLeft = [];
elementsDown = [];
elementsUp = [];


%define side node indices
down_nodes = 1:(p+1)*(q+1);
right_nodes = (p+1):(p+1):(p+1)*(q+1)*(r+1);
left_nodes = 1:(p+1):(1+(p+1)*(q+1)*(r+1)-p);
up_nodes = 1+(p+1)*(q+1)*r:(p+1)*(q+1)*(r+1);

front_nodes = zeros(1,(p+1)*(r+1));
back_nodes = zeros(1, (p+1)*(r+1));

for i=1:r+1
    front_nodes((i-1)*(p+1)+1:i*(p+1)) = (p+1)*(q+1)*(i-1)+1:(p+1)*(q+1)*(i-1)+p+1;
    back_nodes((i-1)*(p+1)+1:i*(p+1)) = i*(p+1)*(q+1)-p:i*(p+1)*(q+1);
end

for elemIndex=1:length(PHTelem)
    if (isempty(PHTelem(elemIndex).children))
        
        %front face
        if (face == 1) && (PHTelem(elemIndex).vertex(2)==0)
            %add the element to the sorted list
            elementsFront = [elementsFront, elemIndex];
            %store the left-lower corner in the verticesFront array
            verticesFront = [verticesFront; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(3)];
        end
        
        %right face
        if (face == 2) && (PHTelem(elemIndex).vertex(4)==1)
            %add the element to the sorted list
            elementsRight = [elementsRight, elemIndex];
            %store the lower-front corner in the verticesRight array
            verticesRight = [verticesRight; PHTelem(elemIndex).vertex(2), PHTelem(elemIndex).vertex(3)];
        end
        
        %back face
        if (face == 3) && (PHTelem(elemIndex).vertex(5)==1)
            %add the element to the sorted list
            elementsBack = [elementsBack, elemIndex];
            %store the left-lower corner in the verticesBack array
            verticesBack = [verticesBack; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(3)];
        end
        
        %left face
        if (face == 4) && (PHTelem(elemIndex).vertex(1)==0)
            %add the element to the sorted list
            elementsLeft = [elementsLeft, elemIndex];
            %store the lower-front corner in the verticesLeft array
            verticesLeft = [verticesLeft; PHTelem(elemIndex).vertex(2), PHTelem(elemIndex).vertex(3)];
        end
        
        %down face
        if (face == 5) && (PHTelem(elemIndex).vertex(3)==0)
            %add the element to the sorted list
            elementsDown = [elementsDown, elemIndex];
            %store the left-front corner in the verticesDown array
            verticesDown = [verticesDown; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(2)];
        end
        
        %up face
        if (face == 6) && (PHTelem(elemIndex).vertex(6)==1)
            %add the element to the sorted list
            elementsUp = [elementsUp, elemIndex];
            %store the lower-front corner in the verticesUp array
            verticesUp = [verticesUp; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(2)];
        end
    end
end

nodes = [];
%sort the elements according to vertex location
switch face
    case 1
        [verticesFront, indexSort] = sortrows(verticesFront);
        %rearrange the elements according to the sorting oder
        elementsFront = elementsFront(indexSort);
        elements = elementsFront;
        for elemIndex = 1:size(verticesFront,1)
            nodes = [nodes, PHTelem(elementsFront(elemIndex)).nodesGlobal(front_nodes)];
        end
    case 2
        [verticesRight, indexSort] = sortrows(verticesRight);
         %rearrange the elements according to the sorting oder
        elementsRight = elementsRight(indexSort);
        elements = elementsRight;
        for elemIndex = 1:size(verticesRight,1)
            nodes = [nodes, PHTelem(elementsRight(elemIndex)).nodesGlobal(right_nodes)];
        end
    case 3
        [verticesBack, indexSort] = sortrows(verticesBack);
         %rearrange the elements according to the sorting oder
        elementsBack = elementsBack(indexSort);
        elements = elementsBack;
        for elemIndex = 1:size(verticesBack,1)
            nodes = [nodes, PHTelem(elementsBack(elemIndex)).nodesGlobal(back_nodes)];
        end
    case 4
        [verticesLeft, indexSort] = sortrows(verticesLeft);
         %rearrange the elements according to the sorting oder
        elementsLeft = elementsLeft(indexSort);
        elements = elementsLeft;
        for elemIndex = 1:size(verticesLeft,1)
            nodes = [nodes, PHTelem(elementsLeft(elemIndex)).nodesGlobal(left_nodes)];
        end
    case 5
        [verticesDown, indexSort] = sortrows(verticesDown);
         %rearrange the elements according to the sorting oder
        elementsDown = elementsDown(indexSort);
        elements = elementsDown;
        for elemIndex = 1:size(verticesDown,1)
            nodes = [nodes, PHTelem(elementsDown(elemIndex)).nodesGlobal(down_nodes)];
        end
    case 6
        [verticesUp, indexSort] = sortrows(verticesUp);
        %rearrange the elements according to the sorting oder
        elementsUp = elementsUp(indexSort);
        elements = elementsUp;
        for elemIndex = 1:size(verticesUp,1)
            nodes = [nodes, PHTelem(elementsUp(elemIndex)).nodesGlobal(up_nodes)];
        end
end


%remove duplicate nodes
%nodes{edgeIndex} = unique_stab(nodes{edgeIndex});
nodes = unique(nodes,'stable');

