function [ elements ] = sortFaceElem( PHTelem, face)
%outputs the elements along the given face, sorted in the order in
%which they appear in the parameter space
%Note: here face is scalar (single value)
%face encoding: 1 - front, 2 - right, 3 - back, 4 - left, 5 - down, 6 - up


%initialized the elements array
elements = [];

%initialize array to store the corners of the edge elements and the edge element indices,
%for sorting purposes
verticesFront = [];
verticesRight = [];
verticesBack = [];
verticesLeft = [];
verticesDown = [];
verticesUp = [];


for elemIndex=1:length(PHTelem)
    if (isempty(PHTelem(elemIndex).children))
        
        %front face
        if (face == 1) && (PHTelem(elemIndex).vertex(2)==0)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the left-lower corner in the verticesFront array
            verticesFront = [verticesFront; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(3)];
        end
        
        %right face
        if (face == 2) && (PHTelem(elemIndex).vertex(4)==1)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the lower-front corner in the verticesRight array
            verticesRight = [verticesRight; PHTelem(elemIndex).vertex(2), PHTelem(elemIndex).vertex(3)];
        end
        
        %back face
        if (face == 3) && (PHTelem(elemIndex).vertex(5)==1)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the left-lower corner in the verticesBack array
            verticesBack = [verticesBack; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(3)];
        end
        
        %left face
        if (face == 4) && (PHTelem(elemIndex).vertex(1)==0)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the lower-front corner in the verticesLeft array
            verticesLeft = [verticesLeft; PHTelem(elemIndex).vertex(2), PHTelem(elemIndex).vertex(3)];
        end
        
        %down face
        if (face == 5) && (PHTelem(elemIndex).vertex(3)==0)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the left-front corner in the verticesDown array
            verticesDown = [verticesDown; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(2)];
        end
        
        %up face
        if (face == 6) && (PHTelem(elemIndex).vertex(6)==1)
            %add the element to the sorted list
            elements = [elements, elemIndex];
            %store the lower-front corner in the verticesUp array
            verticesUp = [verticesUp; PHTelem(elemIndex).vertex(1), PHTelem(elemIndex).vertex(2)];
        end
    end
end

%sort the elements according to vertex location
switch face
    case 1
        [verticesFront, indexSort] = sortrows(verticesFront);
    case 2
        [verticesRight, indexSort] = sortrows(verticesRight);
    case 3
        [verticesBack, indexSort] = sortrows(verticesBack);
    case 4
        [verticesLeft, indexSort] = sortrows(verticesLeft);
    case 5
        [verticesDown, indexSort] = sortrows(verticesDown);
    case 6
        [verticesUp, indexSort] = sortrows(verticesUp);

end

%rearrange the elements according to the sorting oder
elements = elements(indexSort);

end





