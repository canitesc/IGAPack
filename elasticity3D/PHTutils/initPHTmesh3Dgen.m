function [ PHTelem, dimBasis, octupleList ] = initPHTmesh3Dgen( p,q,r, numElemU, numElemV, numElemW )
%initialize the PHT geometry on coarse mesh, with 2*numElemU x 2*numElemV x
%2*numElemW elements (after the first cross insertion)

knotU = [zeros(1,p+1), (1:numElemU)/numElemU, ones(1,p)];
knotV = [zeros(1,q+1), (1:numElemV)/numElemV, ones(1,q)];
knotW = [zeros(1,r+1), (1:numElemW)/numElemW, ones(1,r)];

%repeat the interior knots p-1 times
rep_knotU = linspace(0,1,numElemU+1);
rep_knotU = rep_knotU(2:end-1);
rep_knotU = repmat(rep_knotU,1,p-2);

rep_knotV = linspace(0,1,numElemV+1);
rep_knotV = rep_knotV(2:end-1);
rep_knotV = repmat(rep_knotV,1,q-2);

rep_knotW = linspace(0,1,numElemW+1);
rep_knotW = rep_knotW(2:end-1);
rep_knotW = repmat(rep_knotW,1,r-2);

knotU = sort([knotU, rep_knotU]);
knotV = sort([knotV, rep_knotV]);
knotW = sort([knotW, rep_knotW]);


[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);
[C_w, ~] = bezierExtraction(knotW,r);


%compute the number of basis functions in each direction
lenU = length(knotU)-p-1;
lenV = length(knotV)-q-1;
lenW = length(knotW)-r-1;

dimBasis = lenU*lenV*lenW;
numElements = numElemU*numElemV*numElemW;
octupleList = zeros(numElements, 8);

PHTelem = struct;
%initialize the neighbor connectivity lists
PHTelem.neighbor_left = [];
PHTelem.neighbor_right = [];
PHTelem.neighbor_down = [];
PHTelem.neighbor_up = [];
PHTelem.neighbor_front = [];
PHTelem.neighbor_back = [];
PHTelem.neighbor_up_left = [];
PHTelem.neighbor_down_left = [];
PHTelem.neighbor_up_right = [];
PHTelem.neighbor_down_right = [];
PHTelem.neighbor_up_front = [];
PHTelem.neighbor_down_front = [];
PHTelem.neighbor_up_back = [];
PHTelem.neighbor_down_back = [];
PHTelem.neighbor_left_front = [];
PHTelem.neighbor_right_front = [];
PHTelem.neighbor_left_back = [];
PHTelem.neighbor_right_back = [];


%loop through each element and compute the element-node connectivities
nument = (p+1)*(q+1)*(r+1);
elementCounter = 0;
for k=1:length(knotW)-1
    for j=1:length(knotV)-1
        for i=1:length(knotU)-1
            if (knotU(i+1)>knotU(i)) && (knotV(j+1) >knotV(j)) && (knotW(k+1) > knotW(k))  %the knotspan has non-zero area
                elementCounter = elementCounter + 1;
                PHTelem(elementCounter).parent = [];
                PHTelem(elementCounter).children = [];
                PHTelem(elementCounter).vertex = [knotU(i), knotV(j), knotW(k), knotU(i+1), knotV(j+1), knotW(k+1)];
                
                tcount = 0;
                currow = zeros(1, nument);
                %now we add the nodes from i-p...i in the u
                %direction, j-q...j in the v direction, k-r...k in
                %w direction
                for t3 = k-r:k
                    for t2=j-q:j
                        for t1 = i-p:i
                            tcount = tcount + 1;
                            currow(tcount) = t1+(t2-1)*lenU+(t3-1)*lenU*lenV;
                        end
                    end
                end
                PHTelem(elementCounter).nodes=currow;
                PHTelem(elementCounter).level = 0;
                octupleList(elementCounter,:) = numElements + (8*(elementCounter-1)+1:8*elementCounter);
            end
        end
    end
end


%loop through each element and compute the neighbor lists and Bezier
%extraction operators
indexMatrix = permute(reshape(1:numElements, numElemU, numElemV, numElemW),[2,1,3]);
for k=1:numElemW
    for j=1:numElemV
        for i=1:numElemU
            elementIndex = indexMatrix(j,i,k);
            PHTelem(elementIndex).C = kron(kron(C_w(:,:,k),C_v(:,:,j)),C_u(:,:,i));

            if i>1
                PHTelem(elementIndex).neighbor_left = indexMatrix(j,i-1,k);
                if j>1
                    PHTelem(elementIndex).neighbor_left_front = indexMatrix(j-1,i-1,k);
                end
                if j<numElemV
                    PHTelem(elementIndex).neighbor_left_back = indexMatrix(j+1,i-1,k);
                end
                if k>1
                    PHTelem(elementIndex).neighbor_down_left = indexMatrix(j,i-1,k-1);
                end
                if k<numElemW
                    PHTelem(elementIndex).neighbor_up_left = indexMatrix(j,i-1,k+1);
                end
            end
            if i<numElemU
                PHTelem(elementIndex).neighbor_right = indexMatrix(j,i+1,k);
                if j>1
                    PHTelem(elementIndex).neighbor_right_front = indexMatrix(j-1,i+1,k);
                end
                if j<numElemV
                    PHTelem(elementIndex).neighbor_right_back = indexMatrix(j+1,i+1,k);
                end
                if k>1
                    PHTelem(elementIndex).neighbor_down_right = indexMatrix(j,i+1,k-1);
                end
                if k<numElemW
                    PHTelem(elementIndex).neighbor_up_right = indexMatrix(j,i+1,k+1);
                end
            end
            
            if k>1
                PHTelem(elementIndex).neighbor_down = indexMatrix(j,i,k-1);
                if j>1
                    PHTelem(elementIndex).neighbor_down_front = indexMatrix(j-1,i,k-1);
                end
                if j<numElemV
                    PHTelem(elementIndex).neighbor_down_back = indexMatrix(j+1,i,k-1);
                end
            end
            
            if k<numElemW
                PHTelem(elementIndex).neighbor_up = indexMatrix(j,i,k+1);
                if j>1
                    PHTelem(elementIndex).neighbor_up_front = indexMatrix(j-1,i,k+1);
                end
                if j<numElemV
                    PHTelem(elementIndex).neighbor_up_back = indexMatrix(j+1,i,k+1);
                end
            end
            
            if j>1
                PHTelem(elementIndex).neighbor_front = indexMatrix(j-1,i,k);
            end
            if j<numElemV
                PHTelem(elementIndex).neighbor_back = indexMatrix(j+1,i,k);
            end                        
        end
    end
end

[ PHTelem, dimBasis ] = crossInsert3D( PHTelem, 1:numElements, dimBasis, p, q ,r );


