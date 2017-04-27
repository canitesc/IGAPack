function [ distance ] = makeConformingDist3D( PHTelemA, PHTelemB, GIFTmeshA, GIFTmeshB, elementsA, elementsB, faceA, faceB, p)
%checks that patchA and patchB are conforming (the elements match up) and
%computes distance between the corresponding gauss points on the face
%elements
%Note: the distance should be close to zero (up to round-off errors) if the
%elements match and have correct parametric direction


%make list of coordinates of element vertices along the patch boundaries
numElementsA = length(elementsA);
%numElementsB = length(elementsB);

ngauss_face = p+1;
[gauss_weight_face, gauss_coord_face] = quadrature( ngauss_face, 'GAUSS', 1 );

distance = 0;
%loop over the elements in A and B
for elementIndex = 1:numElementsA
    curElemA = elementsA(elementIndex);
    curElemB = elementsB(elementIndex);
    
    xminA = PHTelemA(curElemA).vertex(1);
    xmaxA = PHTelemA(curElemA).vertex(4);
    yminA = PHTelemA(curElemA).vertex(2);
    ymaxA = PHTelemA(curElemA).vertex(5);
    zminA = PHTelemA(curElemA).vertex(3);
    zmaxA = PHTelemA(curElemA).vertex(6);
    
    xminB = PHTelemB(curElemB).vertex(1);
    xmaxB = PHTelemB(curElemB).vertex(4);
    yminB = PHTelemB(curElemB).vertex(2);
    ymaxB = PHTelemB(curElemB).vertex(5);
    zminB = PHTelemB(curElemB).vertex(3);
    zmaxB = PHTelemB(curElemB).vertex(6);
    
    scalefacA = (xmaxA - xminA)*(ymaxA - yminA)*(zmaxA-zminA)/8;
    %scalefacB = (xmaxB - xminB)*(ymaxB - yminB)*(zmaxB-zminB)/8;
    for jj=1:ngauss_face
        for ii=1:ngauss_face
            
            switch faceA
                case 1 %front
                    u_hatA = gauss_coord_face(ii);
                    v_hatA = -1;
                    w_hatA = gauss_coord_face(jj);
                    
                case 2 %right
                    u_hatA = 1;
                    v_hatA = gauss_coord_face(ii);
                    w_hatA = gauss_coord_face(jj);
                    
                case 3 %back
                    u_hatA = gauss_coord_face(ii);
                    v_hatA = 1;
                    w_hatA = gauss_coord_face(jj);
                    
                case 4 %left
                    u_hatA = -1;
                    v_hatA = gauss_coord_face(ii);
                    w_hatA = gauss_coord_face(jj);
                    
                case 5 %down
                    u_hatA = gauss_coord_face(ii);
                    v_hatA = gauss_coord_face(jj);
                    w_hatA = -1;
                    
                case 6 %up
                    u_hatA = gauss_coord_face(ii);
                    v_hatA = gauss_coord_face(jj);
                    w_hatA = 1;                    
            end
            
            switch faceB
                case 1 %front
                    u_hatB = gauss_coord_face(ii);
                    v_hatB = -1;
                    w_hatB = gauss_coord_face(jj);
                    
                case 2 %right
                    u_hatB = 1;
                    v_hatB = gauss_coord_face(ii);
                    w_hatB = gauss_coord_face(jj);
                    
                case 3 %back
                    u_hatB = gauss_coord_face(ii);
                    v_hatB = 1;
                    w_hatB = gauss_coord_face(jj);
                    
                case 4 %left
                    u_hatB = -1;
                    v_hatB = gauss_coord_face(ii);
                    w_hatB = gauss_coord_face(jj);
                    
                case 5 %down
                    u_hatB = gauss_coord_face(ii);
                    v_hatB = gauss_coord_face(jj);
                    w_hatB = -1;
                    
                case 6 %up
                    u_hatB = gauss_coord_face(ii);
                    v_hatB = gauss_coord_face(jj);
                    w_hatB = 1;                    
            end
           % [u_hatA, v_hatA, w_hatA]
           % [u_hatB, v_hatB, w_hatB]
            
            %evaluate the derivatives of the mapping from parameter
            %space to physical space                                    
            [coordA, dxdxiA] = paramMap3D( GIFTmeshA, u_hatA, v_hatA, w_hatA, xminA, yminA, zminA, xmaxA, ymaxA, zmaxA);
            
            plot3(coordA(1),coordA(2),coordA(3),'.r')
            hold on
            [coordB, ~] = paramMap3D( GIFTmeshB, u_hatB, v_hatB, w_hatB, xminB, yminB, zminB, xmaxB, ymaxB, zmaxB);
            plot3(coordB(1),coordB(2),coordB(3),'.k')
            hold on
            %
%             
%             coordA
%             coordB
%             
%             pause
%             
            J = abs(det(dxdxiA));
            integrand = (coordA(1)-coordB(1))^2 + (coordA(2)-coordB(2))^2 + (coordA(3)-coordB(3))^2;
            distance = distance + integrand * scalefacA * gauss_weight_face(ii).*gauss_weight_face(jj).*J;
        end
    end
    
end
distance = sqrt(distance);
