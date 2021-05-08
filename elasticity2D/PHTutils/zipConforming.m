function [ PHTelem, sizeBasis ] = zipConforming( PHTelem, dimBasis, vertex2patch, edge_list, p, q )
% Connects two conforming patches by changing the nodesGlobal entry using the vertex and
% edge list connectivity
% Input: PHTelem - cell array containing the multi-patch mesh
%        dimBasis - array containing the basis dimension in each patch
%        vertex2patch - cell array containing the vertex-patch connectivity and vertex
%                       corner locations
%        edge_list - structure containing a list of edges (patch interfaces)
%                p - polynomial degree in the u direction (p=3 for "iso" code)
%                q - polynomial degree in the v direction (q=3 for "iso" code)
% Output: PHTelem - updated PHTelem cell array including a nodesGlobal field for each
%                   patch containing the global node numbers
%         sizeBasis - integer scalar representing the dimension of the overall geometry
%                     (related to the size of the global system matrix)


numPatches = length(PHTelem);
corner_indices = [1, p+1, (p+1)*(q+1), (p+1)*q+1];
assignedNodes = cell(1, numPatches);
nodesPattern = cell(1, numPatches);

% initialize nodesPattern in all patches to zero
for patchIndex = 1:numPatches
    nodesPattern{patchIndex} = zeros(1, dimBasis(patchIndex));
end

% Step 1: Loop over the vertices in vertex2patch and assign the global indices as the
% patch corners
sizeBasis = 0;
for vertexindex = 1:length(vertex2patch)
    for i = 1:size(vertex2patch{vertexindex},1)
        patchIndex = vertex2patch{vertexindex}(i,1);
        cornerIndex = vertex2patch{vertexindex}(i,2);
        elementIndex = getElementFromCorner(PHTelem{patchIndex}, cornerIndex);
        localNodes = PHTelem{patchIndex}(elementIndex).nodes(corner_indices(cornerIndex));
        assignedNodes{patchIndex} = [assignedNodes{patchIndex}, localNodes];
        nodesPattern{patchIndex}(localNodes) = vertexindex;
    end
end

% Step 2: Loop over the edges in edge_list and assign the global indices in the interior
% of each edge
sizeBasis = sizeBasis + length(vertex2patch);
edge_fields = fieldnames(edge_list);
for edge_index = 1:length(edge_fields)
    edge_field = edge_fields{edge_index};
    patchA = edge_list.(edge_field)(1);
    edgeA = edge_list.(edge_field)(2);
    nodesAcell = sortEdgeNodesElem( PHTelem{patchA}, edgeA, p, q, 1 );
    nodesA = nodesAcell{1}(2:end-1);
    assignedNodes{patchA} = [assignedNodes{patchA}, nodesA];
    num_new_nodes = length(nodesA);
    newNodeSet = sizeBasis+1:sizeBasis+num_new_nodes;
    nodesPattern{patchA}(nodesA) = newNodeSet;
    % if the edge is in the interior (has length=5)
    if length(edge_list.(edge_field))==5
        %get the nodes on the boundary edge in patchB
        patchB = edge_list.(edge_field)(3);
        edgeB = edge_list.(edge_field)(4);
        flagDir = edge_list.(edge_field)(5);
        nodesBcell = sortEdgeNodesElem( PHTelem{patchB}, edgeB, p, q, flagDir);
        nodesB = nodesBcell{1}(2:end-1);
        assignedNodes{patchB} = [assignedNodes{patchB}, nodesB];
        assert(length(nodesA)==length(nodesB), 'Non-conforming patches encountered.');
        nodesPattern{patchB}(nodesB) = newNodeSet;
    end
    sizeBasis = sizeBasis + length(nodesA);
end

% Step 3: Loop over the patches and assign the interior global indices
for patchIndex = 1:numPatches
   unassignedNodes = setdiff(1:dimBasis(patchIndex), assignedNodes{patchIndex});
   num_new_nodes = length(unassignedNodes);
   newNodeSet = sizeBasis+1:sizeBasis+num_new_nodes;
   nodesPattern{patchIndex}(unassignedNodes) = newNodeSet;
   % loop over the active elements and assign the global node indices
   for iElem = 1:length(PHTelem{patchIndex})
       if isempty(PHTelem{patchIndex}(iElem).children)
           PHTelem{patchIndex}(iElem).nodesGlobal =...
               nodesPattern{patchIndex}(PHTelem{patchIndex}(iElem).nodes);
       end
   end
   sizeBasis = sizeBasis + num_new_nodes;
end


