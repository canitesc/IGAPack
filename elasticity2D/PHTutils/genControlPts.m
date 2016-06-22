function [ controlPts, PHTelem, dimBasis, quadList] = genControlPts( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV )
%generates the controlPts array and the knot vectors for the C1 mesh

dim = 2; %number of physical dimensions

%tolerance for equality tests
toleq = 1e-10;

nrb = nrbmak(coefs,{knotU,knotV});
p_init = nrb.order(1)-1;
q_init = nrb.order(2)-1;
%increase polynomial order
nrb = nrbdegelev(nrb,[p-p_init,q-q_init]);

%TODO: this may need adjustments if the knots are already repeated in the
%initial mesh


%repeat the interior knots p-1 times
rep_knotU = knotU(p+2:end-p-1);
rep_knotU = repmat(rep_knotU,1,p-2);

rep_knotV = knotV(p+2:end-p-1);
rep_knotV = repmat(rep_knotV,1,q-2);

nrb = nrbkntins(nrb,{rep_knotU,rep_knotV});

nrbctrlplot(nrb)
hold on

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
        controlPts(index,:) = [coefs(1,i,j)./coefs(4,i,j), coefs(2,i,j)./coefs(4,i,j), coefs(4,i,j)]; %put the (i,j) node in the coordinate array
    end %for i
end %for j


%do Bezier extraction
numElements = numberElementsU*numberElementsV;

%make elementVertex and elementNode arrays
elementCounter = 0;
PHTelem = struct;
PHTelem.parent = [];
PHTelem.children = [];
PHTelem.level = 0;
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];
quadList = zeros(numElements, 4);
for j=1:length(knotV)-1
    for i=1:length(knotU)-1
        if (abs(knotU(i+1)-knotU(i))>toleq) && (abs(knotV(j+1)-knotV(j))>toleq)
            elementCounter = elementCounter + 1;
            
            PHTelem(elementCounter).vertex = [knotU(i), knotV(j), knotU(i+1), knotV(j+1)];
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
            PHTelem(elementCounter).nodes=currow;
            quadList(elementCounter,:) = numElements + (4*(elementCounter-1)+1:4*elementCounter);
        end
    end
end
[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);
%pause

%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = permute(reshape(1:numElements, numberElementsU, numberElementsV),[2,1]);

for j=1:numberElementsV
    for i=1:numberElementsU
        elementIndex = indexMatrix(j,i);
        PHTelem(elementIndex).C = kron(C_v(:,:,j),C_u(:,:,i));
        
        if i>1
            PHTelem(elementIndex).neighbor_left = indexMatrix(j,i-1);            
        end
        if i<numberElementsU
            PHTelem(elementIndex).neighbor_right = indexMatrix(j,i+1);            
        end
                
        if j>1
            PHTelem(elementIndex).neighbor_down = indexMatrix(j-1,i);
        end
        if j<numberElementsV
            PHTelem(elementIndex).neighbor_up = indexMatrix(j+1,i);
        end
    end
end

%do the initial refinement
[ PHTelem, controlPts, dimBasis ] = crossInsertIso( PHTelem, controlPts, 1:numElements, dimBasis, p, q );
