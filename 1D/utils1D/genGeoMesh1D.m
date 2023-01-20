function [ GeoMesh] = genGeoMesh1D( knotU, coefs, p, numberElementsU )
%generates a GeoMesh object from the knot vectors, control points/weights
%and polynomial degrees given
%uses the same format as NURBS toolbox

dim = 1; %number of physical dimensions

%tolerance for equality tests
toleq = 1e-10;

%the number of control points in the u and v directions
lenU = length(knotU)-p-1;  %number of basis functions in the u direction
numnodes = lenU;      %number of control points
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights

for index=1:lenU % for each node in the x direction
    if size(coefs,1)<4
        coordinates(index,:) = [coefs(1,index), 1]; %put the ith node + weight=1 in the coordinate array
    else
        coordinates(index,:) = [coefs(1,index)./coefs(4,index), coefs(4,index)]; %put the ith node and weight in the coordinate array
    end
end %for i

%do Bezier extraction
C = bezierExtraction(knotU,p);

%make elementVertex and elementNode arrays
elementCounter = 0;
elementVertex = zeros(numberElementsU,2);
elementNode = zeros(numberElementsU, (p+1));


for i=1:length(knotU)-1
    if (abs(knotU(i+1)-knotU(i))>toleq) 
        elementCounter = elementCounter + 1;
        elementVertex(elementCounter, :) = [knotU(i), knotU(i+1)];
        
        %now we add the nodes from i-p...i in the u direction 
        currow = (i-p):i;        
        elementNode(elementCounter,:)=currow;
    end
end


GeoMesh.numberElements = numberElementsU;
GeoMesh.p = p;
GeoMesh.c_net = coordinates;
GeoMesh.C = C;
GeoMesh.elementNode = elementNode;
GeoMesh.elementVertex = elementVertex;

end

