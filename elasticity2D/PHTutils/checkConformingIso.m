function [PHTelem,controlPts,dimBasis,quadList] = checkConformingIso(PHTelem, controlPts, dimBasis, edge_list, p, q, quadList)
% checks that the patches are conforming and if needed makes them conforming
% through mesh refinement using the vertex2patch and edge_list connectivity
% Input: PHTelem - cell array containing the multi-patch mesh
%        controlPts - cell array containing the control points for each patch
%        dimBasis - array containing the basis dimension in each patch
%        edge_list - structure containing a list of edges (patch interfaces)
%                p - polynomial degree in the u direction (p=3 for "iso" code)
%                q - polynomial degree in the v direction (q=3 for "iso" code)
%         quadList - cell array containing the "recovery patches" which are groups of 4
%                    elements of the same refinement level, used in the recovery-based
%                    error estimator in each geometry patch (not to be confused with
%                    the "quads" which represent the initial multipatch geometry from a
%                    quad mesher
% Output: PHTelem - updated PHTelem cell array containing interface-conforming meshes
%                   (after cross-insertions done to ensure interface conformity)
%         controlPts - updated controlPts array for the interface-conforming meshes
%         dimBasis - updated dimBasis array
%         quadList - updated quadList array

keepChecking = true;
while keepChecking
    keepChecking = false;
    edge_fields = fieldnames(edge_list); 
    for edge_index = 1:length(edge_fields)
        edge_field = edge_fields{edge_index};
        % if the edge is in the interior (has length=5)
        if length(edge_list.(edge_field))==5
            %get the elements on the boundary edge in patchA and patchB
            patchA = edge_list.(edge_field)(1);
            patchB = edge_list.(edge_field)(3);
            edgeA = edge_list.(edge_field)(2);
            edgeB = edge_list.(edge_field)(4);
            flagDir = edge_list.(edge_field)(5);
            
            quadListA = quadList{patchA};
            quadListB = quadList{patchB};
            
            [elementsA] = sortEdgeElemIso( PHTelem{patchA}, edgeA, 1);
            [elementsB] = sortEdgeElemIso( PHTelem{patchB}, edgeB, flagDir);
            
            [quadRefA, quadRefB] = makeConformingIso(PHTelem{patchA}, PHTelem{patchB}, elementsA{1},...
                elementsB{1}, edgeA, edgeB, quadListA, quadListB, flagDir);
            indexQuadA = find(quadRefA > 0);
            indexQuadB = find(quadRefB > 0);
            
            if ~isempty(indexQuadA)
                keepChecking = true;
                numNewPatches = length(indexQuadA);
                disp(['In patch ', num2str(patchA), ' refining ',num2str(numNewPatches), ...
                    ' quadruplets to keep conformity with patch ', num2str(patchB)])
                [quadList{patchA}, PHTelem{patchA}, controlPts{patchA}, dimBasis(patchA)] = refineMeshGradedIso(quadRefA, ...
                    quadListA, PHTelem{patchA}, controlPts{patchA}, p, q, dimBasis(patchA));
            end
            
            if ~isempty(indexQuadB)
                keepChecking = true;
                numNewPatches = length(indexQuadB);
                disp(['In patch ', num2str(patchB), ' refining ',num2str(numNewPatches), ...
                    ' quadruplets to keep conformity with patch ', num2str(patchA)])
                [quadList{patchB}, PHTelem{patchB}, controlPts{patchB}, dimBasis(patchB)] = refineMeshGradedIso(quadRefB,...
                    quadListB, PHTelem{patchB}, controlPts{patchB}, p, q, dimBasis(patchB));
            end
        end
    end
end



