function [ PHTmesh, patchList, newC ] = updateMeshPHT( PHTmesh, patchList, indexPatch, patchRef, elementNodes, C, L )
%update the control points and other information for PHTmesh after refinement
%use knot insertion subroutines from NURBS toolbox for calculating the new
%control points/weights

p = PHTmesh.p;
knotVector = PHTmesh.knotVector;
numberElements = PHTmesh.numberElements;
b_net = PHTmesh.b_net;

numNewPatches = length(indexPatch);
numPatches = size(patchList,1);
newKnots = zeros(1, 2*numNewPatches);
elementCounter = 0;
patchCounter = 0;
knotCounter = 0;
newPatchList = zeros(numPatches+numNewPatches,2);

newC = zeros(p+1, p+1, numberElements + 2*numNewPatches);
new_b_net = zeros(size(b_net,1)+2*(p-1)*numNewPatches, 2);

new_b_net(1:2,:) = b_net(1:2, :);
shapeCounter = 2;

for i=1:numPatches
    newPatchList(patchCounter + 1,:) = [elementCounter + 1, elementCounter + 2];
    if patchRef(i)==1 %patch is to be refined
        %add new midpoints to knot list and new elements to patch
        %list
        e = patchList(i,1);
        Ce = C(:,:,e);               
        ni = elementNodes(e,end);                
        xmin = knotVector(ni);
        xmax = knotVector(ni+1);
        xmid = (knotVector(ni+1)+knotVector(ni))/2;
        %calculate the new control points        
        scrt = elementNodes(e,:);
        cpts = b_net(scrt,1);
        w = b_net(scrt, 2);                
        [Ce1, Ce2, newB] = deCasteljau(Ce, xmin, xmax, cpts, w);        
        newC(:,:,elementCounter+1) = Ce1;
        newC(:,:,elementCounter+2) = Ce2;        
        newKnots(knotCounter+1) = xmid;                
        new_b_net(shapeCounter+1:shapeCounter+2*p-4, :) = [newB(:,1), newB(:,2)];
        new_b_net(shapeCounter+2*p-3:shapeCounter+2*p-2, :) = b_net(scrt(end-1:end),:);                        
        
        
        %do the same for second element in the patch
        e = patchList(i,2);
        Ce = C(:,:,e);         
        ni = elementNodes(e,end);
        xmin = knotVector(ni);
        xmax = knotVector(ni+1);
        xmid = (knotVector(ni+1)+knotVector(ni))/2;
        %calculate the new control points   
        scrt = elementNodes(e,:);
        cpts = b_net(scrt,1);
        w = b_net(scrt, 2);
        [Ce1, Ce2, newB] = deCasteljau(Ce, xmin, xmax, cpts, w);        
        
        newC(:,:,elementCounter+3) = Ce1;
        newC(:,:,elementCounter+4) = Ce2;                
        newKnots(knotCounter+2) = xmid;
        newPatchList(patchCounter + 2,:) = [elementCounter + 3, elementCounter + 4];
        
        new_b_net(shapeCounter+2*p-1:shapeCounter+4*p-6, :) = [newB(:,1), newB(:,2)];
        new_b_net(shapeCounter+4*p-5:shapeCounter+4*p-4, :) = b_net(scrt(end-1:end),:);                       
        
        
        %increment counters
        patchCounter = patchCounter + 2;
        elementCounter = elementCounter + 4;
        knotCounter = knotCounter + 2;
        shapeCounter = shapeCounter + 4*(p-1);
    else        
        %add the same local extraction operators to the list
        newC(:,:,elementCounter+1) = C(:,:,patchList(i,1));
        newC(:,:,elementCounter+2) = C(:,:,patchList(i,2));
        
        %add the same control points and weights to the list
        e = patchList(i,1);
        scrt = elementNodes(e,:);
        new_b_net(shapeCounter+1:shapeCounter+p-1,:) = b_net(scrt(3:end), :);
        e = patchList(i,2);
        scrt = elementNodes(e,:);        
        new_b_net(shapeCounter+p:shapeCounter+2*p-2,:) = b_net(scrt(3:end), :);
        
        %increment the counters
        patchCounter = patchCounter + 1;
        elementCounter = elementCounter + 2;
        shapeCounter = shapeCounter + 2*p-2;

    end

end

%repeat the knots
rep_knot = repmat(newKnots,1,p-1);
newKnotVector = sort([knotVector, rep_knot]);

%update the PHTmesh object with the new knots, b_net, numberElements
PHTmesh.knotVector = newKnotVector;
PHTmesh.b_net = new_b_net;
PHTmesh.numberElements = elementCounter;
patchList = newPatchList;
