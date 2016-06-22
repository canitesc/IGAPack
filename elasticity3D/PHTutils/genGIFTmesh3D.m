function [ GIFTmesh ] = genGIFTmesh3D(knotU, knotV, knotW, coefs, p, q, r, numberElementsU, numberElementsV, numberElementsW)
%generates a GIFTmesh object from the knot vectors, control points/weights
%and polynomial degrees given
%uses the same format as NURBS toolbox

dim = 3; %number of physical dimensions
C = zeros((p+1)*(q+1)*(r+1),(p+1)*(q+1)*(r+1),numberElementsU*numberElementsV*numberElementsW);

toleq = 1e-10;

%the number of control points in the u and v directions
lenU = length(knotU)-p-1;  %number of basis functions in the u direction
lenV = length(knotV)-q-1; %number of basis functions in the v direction
lenW = length(knotW)-r-1;
numnodes = lenU*lenV;      %number of control points
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights
index = 0;
for k=1:lenW
    for j=1:lenV %for each node in the y direction
        for i=1:lenU % for each node in the x direction
            index = index + 1; %the index of the coordinate array
            coordinates(index,:) = [coefs(1,i,j,k)./coefs(4,i,j,k), coefs(2,i,j,k)./coefs(4,i,j,k), coefs(3,i,j,k)./coefs(4,i,j,k), coefs(4,i,j,k)]; %put the (i,j,k) node in the coordinate array
        end %for i
    end %for j
end

[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);
[C_w, ~] = bezierExtraction(knotW,r);

%calcuate the Bezier extraction operators
elementCounter = 0;
for k=1:numberElementsW
    for j=1:numberElementsV
        for i=1:numberElementsU
            elementCounter = elementCounter + 1;
            C(:,:,elementCounter) =  kron(kron(C_w(:,:,k),C_v(:,:,j)),C_u(:,:,i));
        end
    end
end
%make elementVertex and elementNode arrays
elementCounter = 0;
elementVertex = zeros(numberElementsU*numberElementsV*numberElementsW,6);
elementNode = zeros(numberElementsU*numberElementsV*numberElementsW, (p+1)*(q+1)*(r+1));

for k=1:length(knotW)-1
    for j=1:length(knotV)-1
        for i=1:length(knotU)-1
            if (abs(knotU(i+1)-knotU(i))>toleq) && (abs(knotV(j+1)-knotV(j))>toleq) && (abs(knotW(k+1)-knotW(k))>toleq)
                elementCounter = elementCounter + 1;
                
                elementVertex(elementCounter, :) = [knotU(i), knotV(j), knotW(k), knotU(i+1), knotV(j+1), knotW(k+1)];
                
                tcount = 0;
                currow = zeros(1, (p+1)*(q+1)*(r+1));
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
                elementNode(elementCounter,:)=currow;
            end
        end
    end
end

GIFTmesh.numberElements = numberElementsU*numberElementsV*numberElementsW;
GIFTmesh.numberElementsU = numberElementsU;
GIFTmesh.numberElementsV = numberElementsV;
GIFTmesh.numberElementsW = numberElementsW;

GIFTmesh.p = p;
GIFTmesh.q = q;
GIFTmesh.r = r;
GIFTmesh.c_net = coordinates;
GIFTmesh.C = C;
GIFTmesh.elementNode = elementNode;
GIFTmesh.elementVertex =elementVertex;

end

