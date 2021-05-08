function [vertices, vertex2patch, patch2vertex] = genVertex2Patch2D(patch_list, controlPts, p, q)
%  Generates a list of vertex corners and the connectivity between the corner
%     vertices and the patches
% 
%     Parameters
%     ----------
%     patch_list : cell array of PHTelem structs
%     controlPts :  cell array containing the control points in each patch
%     p : polynomial degree in u direction
%     q : polynomial degree in v direction
% 
%     Returns
%     -------
%     vertices : (array of length N) list of corner vertices containing the x, y and z 
%                coordinate
%     vertex2patch : (list of arrays of length N) list where each entry is an
%                     array, with the first column the patch indices that contain
%                     each vertex and in the second column the location of the
%                     corner of each vertex in a patch with the encoding 
%                          1 - down-left  (u=0, v=0), 
%                          2 - down-right (u=1, v=0), 
%                          3 - up-right   (u=1, v=1), 
%                          4 - up-left    (u=0, v=1)
%                          
%     patch2vertex : (2D arrays with N rows and 4 columns) array which gives the patch to 
%                     vertex connectivity, where the vertices are in counter-
%                     clockwise order starting from the origin of the parameter 
%                     space (u=0, v=0)

numPatches = length(patch_list);
vertices = [];
vertex2patch = {};
patch2vertex = cell(1, numPatches);

tol_eq = 1e-10;
% loop over all the patches and elements and add the corner vertices
for i=1:numPatches
    cur_patch = patch_list{i};
    patch_entry = zeros(1, 4);
    num_elem = length(cur_patch);
    cpts = controlPts{i};
    for j = 1:num_elem
        if isempty(cur_patch(j).children)
            elem_vertex = cur_patch(j).vertex;
            elem_node = cur_patch(j).nodes;
            if abs(elem_vertex(1))<tol_eq && abs(elem_vertex(2))<tol_eq
                % get the vertex coordinates for the down_left element 
                encoding = 1;
                node_index = getCornerNode2D([p,q], 'down_left');
                
                [vertices, vertex2patch, patch_entry] = updateIndices(elem_node, cpts, ...
                    vertices, vertex2patch, encoding, patch_entry, node_index, i);
            end
            if abs(elem_vertex(3)-1)<tol_eq && abs(elem_vertex(2))<tol_eq
                % get the vertex coordinates for the down_right element 
                encoding = 2;
                node_index = getCornerNode2D([p,q], 'down_right');
                [vertices, vertex2patch, patch_entry] = updateIndices(elem_node, cpts, ...
                    vertices, vertex2patch, encoding, patch_entry, node_index, i);
            end
            if abs(elem_vertex(3)-1)<tol_eq && abs(elem_vertex(4)-1)<tol_eq
                % get the vertex coordinates for the up_right element 
                encoding = 3;
                node_index = getCornerNode2D([p,q], 'up_right');
                [vertices, vertex2patch, patch_entry] = updateIndices(elem_node, cpts, ...
                    vertices, vertex2patch, encoding, patch_entry, node_index, i);
            end
            if abs(elem_vertex(1))<tol_eq && abs(elem_vertex(4)-1)<tol_eq
                % get the vertex coordinates for the up_left element 
                encoding = 4;
                node_index = getCornerNode2D([p,q], 'up_left');
                [vertices, vertex2patch, patch_entry] = updateIndices(elem_node, cpts, ...
                    vertices, vertex2patch, encoding, patch_entry, node_index, i);
            end
        end
    end
    patch2vertex{i} = patch_entry;
    
end

end

function [vertices, vertex2patch, patch_entry] = updateIndices(elem_node, ...
                        cpts, vertices, vertex2patch, encoding,...
                        patch_entry, node_index, i)
cpt_index = elem_node(node_index);
vertex = cpts(cpt_index, :);
vert_index = findVertex(vertices, vertex);
if isnan(vert_index)
    % create new entry in the vertices and vertex2patch list
   vertices = [vertices; vertex];
   row_entry = [i, encoding];
   vertex2patch{end+1} = row_entry;
   vert_index = size(vertices,1);
else
    % update the vertex2patch list with the current patch
    row_entry = [i, encoding];
    vertex2patch{vert_index} = [vertex2patch{vert_index}; row_entry];    
end
patch_entry(encoding) = vert_index;

   
end
    
