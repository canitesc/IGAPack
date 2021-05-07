function  checkPHTmesh( PHTelem, p, q )
%checks whether the nodes PHTelem have extra nodes

for indexPatch =1:length(PHTelem)
    for indexElem = 1:length(PHTelem{indexPatch})
        if length(PHTelem{indexPatch}(indexElem).nodes) ~= (p+1)*(q+1)
            disp('Extra nodes encountered!')
            indexPatch
            indexElem
            PHTelem{indexPatch}(indexElem).nodes
        end
    end

end

