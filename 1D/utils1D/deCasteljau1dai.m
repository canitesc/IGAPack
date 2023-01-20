function [ Ce1, Ce2, in1, in2, dimBasis ] = deCasteljau1dai( Ce, knotU1, knotU2, p, elmIn, dimBasis )
% Split element with Bezier extraction operator Ce into 2 elements with
% Bezier extraction operator Ce1, Ce2
%calculates element_nod indices for new elements


alpha = floor((p-1)/2);

num_rows = size(Ce,1);
num_cols = size(Ce,2);

Ce1 = zeros(num_rows, num_cols);
Ce2 = zeros(num_rows, num_cols);

knotUmid = (knotU1+knotU2)/2;
newKnotU = [zeros(1,p+1), knotU1*ones(1,p-alpha), knotUmid*ones(1,p-alpha), knotU2*ones(1,p-alpha), ones(1, p+1)];
[newC_u, ~] = bezierExtraction(newKnotU, p);

%calculate the tensor product Bezier extraction operators on the two
%new elements
newC1 = newC_u(:,:,2);
newC2 = newC_u(:,:,3);

%copy the basis function indices to the children elements
in1 = elmIn;
in2 = elmIn;

%do the deCasteljau algorithm for each row
for i=1:size(Ce,1)
    cur_row = Ce(i,:);
        
    [temp1, temp2] = deCasteljau1d(cur_row);
    temp = [temp1, temp2];
    
    %zero out entries coresponding to new vertices
    temp((alpha+2):2*p-alpha+1) = 0;
    
    Ce1(i,:) = temp(1:p+1);
    Ce2(i,:) = temp(p+2:end);
end

%add in the new basis functions
[ west, center, east ] = getCornerIndices1D( p );
[newBasisSet, dimBasis ] = newBasisIndices1D( elmIn, p, dimBasis );

Ce1(center,:) = newC1(center,:);
Ce1(east,:) = newC1(east,:);

Ce2(west,:) = newC2(west,:);
Ce2(center,:) = newC2(center,:);


%over-write the element node indices
in1 = newBasisSet{1};
in2 = newBasisSet{2};
