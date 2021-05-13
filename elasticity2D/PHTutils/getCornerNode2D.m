function index = getCornerNode2D(deg, position)
% Gets the corner node position for an element of given degree
% 
%     Parameters
%     ----------
%     deg : (array like of length 2) 
%         degree in the u and v direction (p and q)
%     position : (string)
%         position of the vertex can be either 
%             "down_left"  (u=0, v=0), 
%             "down_right" (u=1, v=0), 
%             "up_right"   (u=1, v=1), 
%             "up_left"    (u=0, v=1)
% 
%     Returns
%     -------
%     index : (int) index of the node in the elem_node array of length (p+1)*(q+1)

if strcmp(position, 'down_left')
    index=1;
elseif strcmp(position, 'down_right')
    index = deg(1)+1;
elseif strcmp(position, 'up_right')
    index = (deg(1)+1)*(deg(2)+1);
elseif strcmp(position, 'up_left')
    index = (deg(1)+1)*deg(2)+1;
else
    error('Wrong position given')
end


