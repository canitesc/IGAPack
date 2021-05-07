function [edge_list] = genEdgeList(patch2vertex)
% Generates a list of patch boundaries (edges) 
% 
%     Parameters
%     ----------    
%     patch2vertex: (list of arrays) patch to vertex connectivity
% 
%     Returns
%     -------
%     edge_list : dict where each entry contains either an array
%                 of length 5 or an array of length 2. If it is an array of length
%                 5, it is of the form [patchA, sideA, patchB, sideB, direction_flag]
%                 where patchA and patchB are the patch indices sharing the edge, 
%                 direction_flag = 1 if the parameteric direction matches along the
%                 edge and direction_flag = -1 if the parametric direction is
%                 reversed, sideA, sideB give the orientation of the side on patchA
%                 and patchB using the encoding 
%                   1 - down (v=0), 
%                   2 - right (u=1),
%                   3 - up (v=1),
%                   4 - left (u=0). 
%                 If it is an array of length 2, it is of
%                 the form [patch, side], which shows the patch index (for a boundary
%                 patch) and the side orientation as before. The key name is of
%                 the form "n_index1_n_index2" where index1 and index2 are the
%                 indices of the endpoint vertices, and index1<index2. 

edge_list = struct;
num_patches = length(patch2vertex);

for i=1:num_patches
    for j=1:length(patch2vertex{i})
        if j<length(patch2vertex{i})
            pt_a = min(patch2vertex{i}(j), patch2vertex{i}(j+1));
            pt_b = max(patch2vertex{i}(j), patch2vertex{i}(j+1));
        else
            pt_a = min(patch2vertex{i}(j), patch2vertex{i}(1));
            pt_b = max(patch2vertex{i}(j), patch2vertex{i}(1));
        end
        edge_field = strcat('n_', num2str(pt_a), '_n_', num2str(pt_b));
        if isfield(edge_list, edge_field)
            edge_list.(edge_field)(3) = i;
            edge_list.(edge_field)(4) = j;
            % rearrange the order of the quad vertices so that the order
            % matches that of the parametric direction on each edge
            temp_quad_a = patch2vertex{edge_list.(edge_field)(1)}([1,2,4,3]);
            temp_quad_b = patch2vertex{i}([1,2,4,3]);
            [~, indx_pt_a_quad_a] = ismember(pt_a, temp_quad_a);
            [~, indx_pt_a_quad_b] = ismember(pt_a, temp_quad_b);
            [~, indx_pt_b_quad_a] = ismember(pt_b, temp_quad_a);
            [~, indx_pt_b_quad_b] = ismember(pt_b, temp_quad_b);
            % if the point locations on each quad are both in increasing
            % order or both in decreasing order, then the parameteric
            % direction is matching, otherwise it is reversed
            if (indx_pt_a_quad_a > indx_pt_b_quad_a && indx_pt_a_quad_b > indx_pt_b_quad_b) ||...
                (indx_pt_a_quad_a < indx_pt_b_quad_a && indx_pt_a_quad_b < indx_pt_b_quad_b)
                edge_list.(edge_field)(5) = 1;
            else
                edge_list.(edge_field)(5) = -1;
            end
        else
            edge_list.(edge_field)=[i,j];
        end
    end
end


end

