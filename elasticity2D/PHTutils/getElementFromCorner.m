function [ element ] = getElementFromCorner( PHTelem, cornerIndex)
% Outputs the element corresponding to the corner index in the given PHTelem patch
% Input: PHTelem - the PHT elem geometry for a single patch
%        cornerIndex - the given corner index using the encoding 
%                   1 - down-left  (u=0, v=0), 
%                   2 - down-right (u=1, v=0), 
%                   3 - up-right   (u=1, v=1), 
%                   4 - up-left    (u=0, v=1)
% Output: The element index which is active (has no children) corresponding to the given
%           corner index

for elemIndex=1:length(PHTelem)
    if isempty(PHTelem(elemIndex).children)
        switch cornerIndex
            case 1
                if isempty(PHTelem(elemIndex).neighbor_down) && isempty(PHTelem(elemIndex).neighbor_left)
                    element = elemIndex;
                    break
                end
            case 2
                if isempty(PHTelem(elemIndex).neighbor_down) && isempty(PHTelem(elemIndex).neighbor_right)
                    element = elemIndex;
                    break
                end
            case 3
                if isempty(PHTelem(elemIndex).neighbor_up) && isempty(PHTelem(elemIndex).neighbor_right)
                    element = elemIndex;
                    break
                end
            case 4
                if isempty(PHTelem(elemIndex).neighbor_up) && isempty(PHTelem(elemIndex).neighbor_left)
                    element = elemIndex;
                    break
                end
            otherwise
                error('Wrong corner index given!')
        end
        
    end
end
