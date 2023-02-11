function [ controlPts, IGAelem, dimBasis, nrb] = genControlPtsNURBS( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV )
%generates the controlPts array and the knot vectors for the NURBS mesh
% no initial repeating of control points

dim = 2; %number of physical dimensions

%tolerance for equality tests
toleq = 1e-10;

nrb = geometry_deg_elevate(knotU, knotV, coefs, p, q);
knotU = nrb.knots{1};
knotV = nrb.knots{2};
coefs = nrb.coefs;

%the number of control points in the u and v directions
lenU = length(knotU)-p-1;  %number of basis functions in the u direction
lenV = length(knotV)-q-1;  %number of basis functions in the v direction
dimBasis = lenU*lenV;      %number of control points
controlPts = zeros(dimBasis, dim+1);  %allocate one more column for the weights

for j=1:lenV %for each node in the y direction
    for i=1:lenU % for each node in the x direction
        index = (j-1)*lenU + i; %the index of the coordinate array
        controlPts(index,:) = [coefs(1,i,j)./coefs(4,i,j), coefs(2,i,j)./coefs(4,i,j),...
            coefs(4,i,j)]; %put the (i,j) node in the coordinate array
    end %for i
end %for j


%do Bezier extraction
numElements = numberElementsU*numberElementsV;

%make elementVertex and elementNode arrays
elementCounter = 0;
IGAelem = struct;

IGAelem.neighbor_left = [];
IGAelem.neighbor_right = [];
IGAelem.neighbor_down = [];
IGAelem.neighbor_up = [];
for j=1:length(knotV)-1
    for i=1:length(knotU)-1
        if (abs(knotU(i+1)-knotU(i))>toleq) && (abs(knotV(j+1)-knotV(j))>toleq)
            elementCounter = elementCounter + 1;
            IGAelem(elementCounter).parent = [];
            IGAelem(elementCounter).children = [];
            IGAelem(elementCounter).level = 0;
            IGAelem(elementCounter).vertex = [knotU(i), knotV(j), knotU(i+1), knotV(j+1)];
            tcount = 0;
            currow = zeros(1, (p+1)*(q+1));
            %now we add the nodes from i-p...i in the u direction and
            %j-q...j in the v direction
            for t2=j-q:j
                for t1 = i-p:i
                    tcount = tcount + 1;
                    currow(tcount) = t1+(t2-1)*lenU;
                end
            end
            IGAelem(elementCounter).nodes=currow;
        end
    end
end
[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);


%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = permute(reshape(1:numElements, numberElementsU, numberElementsV),[2,1]);

for j=1:numberElementsV
    for i=1:numberElementsU
        elementIndex = indexMatrix(j,i);
        IGAelem(elementIndex).C = kron(C_v(:,:,j),C_u(:,:,i));
        
        if i>1
            IGAelem(elementIndex).neighbor_left = indexMatrix(j,i-1);            
        end
        if i<numberElementsU
            IGAelem(elementIndex).neighbor_right = indexMatrix(j,i+1);            
        end
                
        if j>1
            IGAelem(elementIndex).neighbor_down = indexMatrix(j-1,i);
        end
        if j<numberElementsV
            IGAelem(elementIndex).neighbor_up = indexMatrix(j+1,i);
        end
    end
end
