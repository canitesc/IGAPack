function [ bez, nodes, dimBasis ] = deCasteljau3d( Ce, knotU1, knotU2, knotV1, knotV2, knotW1, knotW2, p, q, r, newVert, knotUl, knotUr, knotVd, knotVu, knotWb, knotWt, elmIn, dimBasis )
% Split element with Bezier extraction operator Ce into 8 elements with
% Bezier extraction operator bez.Ce1, bez.Ce2, ..., bez.Ce8
%newVert: new Vertices
%NO automatic detection of extra basis vertices (saved entries)!
%calculates element_nod indices for new elements
%encoding scheme newVert: 1-front, 2-right, 3-back, 4-left, 5-down, 6-up, 7-center  
%8-up_left, 9-down_left, 


num_rows = size(Ce,1);
num_cols = size(Ce,2);

bez.Ce1 = zeros(num_rows, num_cols);
bez.Ce2 = zeros(num_rows, num_cols);
bez.Ce3 = zeros(num_rows, num_cols);
bez.Ce4 = zeros(num_rows, num_cols);
bez.Ce5 = zeros(num_rows, num_cols);
bez.Ce6 = zeros(num_rows, num_cols);
bez.Ce7 = zeros(num_rows, num_cols);
bez.Ce8 = zeros(num_rows, num_cols);

knotUmid = (knotU1+knotU2)/2;
knotVmid = (knotV1+knotV2)/2;
knotWmid = (knotW1+knotW2)/2;

newKnotU = [knotUl*ones(1,p+1), knotU1*ones(1,p-1), knotUmid*ones(1,p-1), knotU2*ones(1,p-1), knotUr*ones(1, p+1)];
newKnotV = [knotVd*ones(1,q+1), knotV1*ones(1,q-1), knotVmid*ones(1,q-1), knotV2*ones(1,q-1), knotVu*ones(1, q+1)];
newKnotW = [knotWb*ones(1,r+1), knotW1*ones(1,r-1), knotWmid*ones(1,r-1), knotW2*ones(1,r-1), knotWt*ones(1, q+1)];

[newC_u, ~] = bezierExtraction(newKnotU, p);
[newC_v, ~] = bezierExtraction(newKnotV, q);
[newC_w, ~] = bezierExtraction(newKnotW, r);

%calculate the tensor product Bezier extraction operators on the eight
%new elements
newC1 = kron(kron(newC_w(:,:,2),newC_v(:,:,2)),newC_u(:,:,2));
newC2 = kron(kron(newC_w(:,:,2),newC_v(:,:,2)),newC_u(:,:,3));
newC3 = kron(kron(newC_w(:,:,2),newC_v(:,:,3)),newC_u(:,:,2));
newC4 = kron(kron(newC_w(:,:,2),newC_v(:,:,3)),newC_u(:,:,3));
newC5 = kron(kron(newC_w(:,:,3),newC_v(:,:,2)),newC_u(:,:,2));
newC6 = kron(kron(newC_w(:,:,3),newC_v(:,:,2)),newC_u(:,:,3));
newC7 = kron(kron(newC_w(:,:,3),newC_v(:,:,3)),newC_u(:,:,2));
newC8 = kron(kron(newC_w(:,:,3),newC_v(:,:,3)),newC_u(:,:,3));

%copy the basis function indices to the children elements
nodes.in1 = elmIn;
nodes.in2 = elmIn;
nodes.in3 = elmIn;
nodes.in4 = elmIn;
nodes.in5 = elmIn;
nodes.in6 = elmIn;
nodes.in7 = elmIn;
nodes.in8 = elmIn;

%do tensor product deCasteljau
for i=1:size(Ce,1)
    temp = zeros(2*(q+1), 2*(p+1), 2*(r+1));
    cur_row = Ce(i,:);
    
    cur_row_cube = permute(reshape(cur_row,p+1,q+1,r+1),[2,1,3]);       
    %do the 1d deCasteljau algorithm in the row direction for each slice
    for j=1:q+1
        for k=1:r+1
            [temp1, temp2] = deCasteljau1d(cur_row_cube(j,:,k));
            temp(j,:,k) = [temp1, temp2];
        end
    end
    %do the 1d deCasteljau algorithm in the column direction for each slice
    for j=1:2*(p+1)
        for k=1:r+1
            [temp1, temp2] = deCasteljau1d(temp(1:(q+1),j,k)');
            temp(:,j,k) = [temp1, temp2];
        end
    end    
    
    %do the 1d deCasteljau algorithm in the slice direction for each row
    %and column
    for j=1:2*(q+1)
        for k=1:2*(p+1)
            [temp1, temp2] = deCasteljau1d(squeeze(temp(j,k,1:(r+1)))');         
            temp(j,k,:) = [temp1, temp2];
        end
    end                    

    %zero out entries coresponding to new basis vertices
    for j=1:length(newVert)
        switch newVert(j)
            case 1
                temp(1:q-1,3:2*p,3:2*r) = 0;
            case 2
                temp(3:2*q,end-p+2:end,3:2*r) = 0;
            case 3
                temp(end-q+2:end,3:2*p,3:2*r) = 0;
            case 4
                temp(3:2*q,1:p-1,3:2*r) = 0;
            case 5
                temp(3:2*q,3:2*p,1:r-1) = 0;
            case 6
                temp(3:2*q,3:2*p,end-r+2:end) = 0;
            case 7
                temp(3:2*q,3:2*p,3:2*r) = 0;
            case 8
                temp(3:2*q,1:p-1,end-r+2:end) = 0;
            case 9
                temp(3:2*q,1:p-1,1:r-1) = 0;
            case 10
                temp(3:2*q,end-p+2:end,end-r+2:end) = 0;
            case 11
                temp(3:2*q,end-p+2:end,1:r-2) = 0;
            case 12
                temp(1:q-1,3:2*p,end-r+2:end) = 0;
            case 13
                temp(1:q-1,3:2*p,1:r-1) = 0;
            case 14
                temp(end-q+2:end,3:2*p,end-r+2:end) = 0;
            case 15
                temp(end-q+2:end,3:2*p,1:r-1) = 0;
            case 16
                temp(1:q-1,1:p-1,3:2*r) = 0;
            case 17
                temp(1:q-1,end-p+2:end,3:2*r) = 0;
            case 18
                temp(end-q+2:end,1:p-1,3:2*r) = 0;
            case 19
                temp(end-q+2:end,end-p+2:end,3:2*r) = 0;                
        end
    end    
    bez.Ce1(i,:) = reshape(permute(temp(1:q+1,1:p+1,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce2(i,:) = reshape(permute(temp(1:q+1,p+2:end,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce3(i,:) = reshape(permute(temp(q+2:end,1:p+1,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce4(i,:) = reshape(permute(temp(q+2:end,p+2:end,1:r+1),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce5(i,:) = reshape(permute(temp(1:q+1,1:p+1,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce6(i,:) = reshape(permute(temp(1:q+1,p+2:end,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce7(i,:) = reshape(permute(temp(q+2:end,1:p+1,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
    bez.Ce8(i,:) = reshape(permute(temp(q+2:end,p+2:end,r+2:end),[2,1,3]),1,(p+1)*(q+1)*(r+1));
                          
end

%add in the new basis functions
%define corners indices
reg = getRegionIndices3D( p,q,r );

[newBasisSet, dimBasis ] = newBasisIndices3D( newVert, elmIn, p, q, r, dimBasis, reg );
%dimBasis

for j=1:length(newVert)
    switch newVert(j)
        case 1           
            bez.Ce1(reg.se_mid,:) = newC1(reg.se_mid,:);
            bez.Ce1(reg.south_mid,:) = newC1(reg.south_mid,:);
            bez.Ce2(reg.south_mid,:) = newC2(reg.south_mid,:);
            bez.Ce2(reg.sw_mid,:) = newC2(reg.sw_mid,:);        
            bez.Ce1(reg.se_top,:) = newC1(reg.se_top,:);
            bez.Ce1(reg.south_top,:) = newC1(reg.south_top,:);
            bez.Ce2(reg.south_top,:) = newC2(reg.south_top,:);
            bez.Ce2(reg.sw_top,:) = newC2(reg.sw_top,:);
            bez.Ce5(reg.se_low,:) = newC5(reg.se_low,:);
            bez.Ce5(reg.south_low,:) = newC5(reg.south_low,:);
            bez.Ce6(reg.south_low,:) = newC6(reg.south_low,:);
            bez.Ce6(reg.sw_low,:) = newC6(reg.sw_low,:);                        
            bez.Ce5(reg.se_mid,:) = newC5(reg.se_mid,:);
            bez.Ce5(reg.south_mid,:) = newC5(reg.south_mid,:);
            bez.Ce6(reg.south_mid,:) = newC6(reg.south_mid,:);
            bez.Ce6(reg.sw_mid,:) = newC6(reg.sw_mid,:);                        
        case 2
            bez.Ce2(reg.ne_mid,:) = newC2(reg.ne_mid,:);
            bez.Ce2(reg.east_mid,:) = newC2(reg.east_mid,:);
            bez.Ce4(reg.east_mid,:) = newC4(reg.east_mid,:);
            bez.Ce4(reg.se_mid,:) = newC4(reg.se_mid,:);        
            bez.Ce2(reg.ne_top,:) = newC2(reg.ne_top,:);
            bez.Ce2(reg.east_top,:) = newC2(reg.east_top,:);
            bez.Ce4(reg.east_top,:) = newC4(reg.east_top,:);
            bez.Ce4(reg.se_top,:) = newC4(reg.se_top,:);
            bez.Ce6(reg.ne_low,:) = newC6(reg.ne_low,:);
            bez.Ce6(reg.east_low,:) = newC6(reg.east_low,:);
            bez.Ce8(reg.east_low,:) = newC8(reg.east_low,:);
            bez.Ce8(reg.se_low,:) = newC8(reg.se_low,:);                     
            bez.Ce6(reg.ne_mid,:) = newC6(reg.ne_mid,:);
            bez.Ce6(reg.east_mid,:) = newC6(reg.east_mid,:);
            bez.Ce8(reg.east_mid,:) = newC8(reg.east_mid,:);
            bez.Ce8(reg.se_mid,:) = newC8(reg.se_mid,:);                     
        case 3
            bez.Ce3(reg.ne_mid,:) = newC3(reg.ne_mid,:);
            bez.Ce3(reg.north_mid,:) = newC3(reg.north_mid,:);
            bez.Ce4(reg.north_mid,:) = newC4(reg.north_mid,:);
            bez.Ce4(reg.nw_mid,:) = newC4(reg.nw_mid,:);        
            bez.Ce3(reg.ne_top,:) = newC3(reg.ne_top,:);
            bez.Ce3(reg.north_top,:) = newC3(reg.north_top,:);
            bez.Ce4(reg.north_top,:) = newC4(reg.north_top,:);
            bez.Ce4(reg.nw_top,:) = newC4(reg.nw_top,:);
            bez.Ce7(reg.ne_low,:) = newC7(reg.ne_low,:);
            bez.Ce7(reg.north_low,:) = newC7(reg.north_low,:);
            bez.Ce8(reg.north_low,:) = newC8(reg.north_low,:);
            bez.Ce8(reg.nw_low,:) = newC8(reg.nw_low,:);                     
            bez.Ce7(reg.ne_mid,:) = newC7(reg.ne_mid,:);
            bez.Ce7(reg.north_mid,:) = newC7(reg.north_mid,:);
            bez.Ce8(reg.north_mid,:) = newC8(reg.north_mid,:);
            bez.Ce8(reg.nw_mid,:) = newC8(reg.nw_mid,:);                            
        case 4
            bez.Ce1(reg.nw_mid,:) = newC1(reg.nw_mid,:);
            bez.Ce1(reg.west_mid,:) = newC1(reg.west_mid,:);
            bez.Ce3(reg.west_mid,:) = newC3(reg.west_mid,:);
            bez.Ce3(reg.sw_mid,:) = newC3(reg.sw_mid,:);        
            bez.Ce1(reg.nw_top,:) = newC1(reg.nw_top,:);
            bez.Ce1(reg.west_top,:) = newC1(reg.west_top,:);
            bez.Ce3(reg.west_top,:) = newC3(reg.west_top,:);
            bez.Ce3(reg.sw_top,:) = newC3(reg.sw_top,:);
            bez.Ce5(reg.nw_low,:) = newC5(reg.nw_low,:);
            bez.Ce5(reg.west_low,:) = newC5(reg.west_low,:);
            bez.Ce7(reg.west_low,:) = newC7(reg.west_low,:);
            bez.Ce7(reg.sw_low,:) = newC7(reg.sw_low,:);                     
            bez.Ce5(reg.nw_mid,:) = newC5(reg.nw_mid,:);
            bez.Ce5(reg.west_mid,:) = newC5(reg.west_mid,:);
            bez.Ce7(reg.west_mid,:) = newC7(reg.west_mid,:);
            bez.Ce7(reg.sw_mid,:) = newC7(reg.sw_mid,:);                           
        case 5
            bez.Ce1(reg.center_low,:) = newC1(reg.center_low,:);
            bez.Ce1(reg.east_low,:) = newC1(reg.east_low,:);
            bez.Ce1(reg.north_low,:) = newC1(reg.north_low,:);
            bez.Ce1(reg.ne_low,:) = newC1(reg.ne_low,:);        
            bez.Ce2(reg.center_low,:) = newC2(reg.center_low,:);
            bez.Ce2(reg.west_low,:) = newC2(reg.west_low,:);
            bez.Ce2(reg.north_low,:) = newC2(reg.north_low,:);
            bez.Ce2(reg.nw_low,:) = newC2(reg.nw_low,:);
            bez.Ce3(reg.center_low,:) = newC3(reg.center_low,:);
            bez.Ce3(reg.south_low,:) = newC3(reg.south_low,:);
            bez.Ce3(reg.east_low,:) = newC3(reg.east_low,:);
            bez.Ce3(reg.se_low,:) = newC3(reg.se_low,:);                     
            bez.Ce4(reg.center_low,:) = newC4(reg.center_low,:);
            bez.Ce4(reg.south_low,:) = newC4(reg.south_low,:);
            bez.Ce4(reg.west_low,:) = newC4(reg.west_low,:);
            bez.Ce4(reg.sw_low,:) = newC4(reg.sw_low,:);               
        case 6
            bez.Ce5(reg.center_top,:) = newC5(reg.center_top,:);
            bez.Ce5(reg.east_top,:) = newC5(reg.east_top,:);
            bez.Ce5(reg.north_top,:) = newC5(reg.north_top,:);
            bez.Ce5(reg.ne_top,:) = newC5(reg.ne_top,:);        
            bez.Ce6(reg.center_top,:) = newC6(reg.center_top,:);
            bez.Ce6(reg.west_top,:) = newC6(reg.west_top,:);
            bez.Ce6(reg.north_top,:) = newC6(reg.north_top,:);
            bez.Ce6(reg.nw_top,:) = newC6(reg.nw_top,:);
            bez.Ce7(reg.center_top,:) = newC7(reg.center_top,:);
            bez.Ce7(reg.south_top,:) = newC7(reg.south_top,:);
            bez.Ce7(reg.east_top,:) = newC7(reg.east_top,:);
            bez.Ce7(reg.se_top,:) = newC7(reg.se_top,:);                     
            bez.Ce8(reg.center_top,:) = newC8(reg.center_top,:);
            bez.Ce8(reg.south_top,:) = newC8(reg.south_top,:);
            bez.Ce8(reg.west_top,:) = newC8(reg.west_top,:);
            bez.Ce8(reg.sw_top,:) = newC8(reg.sw_top,:);     
        case 7
            %mid_slice in elements 1,2,3,4
            bez.Ce1(reg.center_mid,:) = newC1(reg.center_mid,:);
            bez.Ce1(reg.east_mid,:) = newC1(reg.east_mid,:);
            bez.Ce1(reg.north_mid,:) = newC1(reg.north_mid,:);
            bez.Ce1(reg.ne_mid,:) = newC1(reg.ne_mid,:);        
            bez.Ce2(reg.center_mid,:) = newC2(reg.center_mid,:);
            bez.Ce2(reg.west_mid,:) = newC2(reg.west_mid,:);
            bez.Ce2(reg.north_mid,:) = newC2(reg.north_mid,:);
            bez.Ce2(reg.nw_mid,:) = newC2(reg.nw_mid,:);
            bez.Ce3(reg.center_mid,:) = newC3(reg.center_mid,:);
            bez.Ce3(reg.south_mid,:) = newC3(reg.south_mid,:);
            bez.Ce3(reg.east_mid,:) = newC3(reg.east_mid,:);
            bez.Ce3(reg.se_mid,:) = newC3(reg.se_mid,:);                     
            bez.Ce4(reg.center_mid,:) = newC4(reg.center_mid,:);
            bez.Ce4(reg.south_mid,:) = newC4(reg.south_mid,:);
            bez.Ce4(reg.west_mid,:) = newC4(reg.west_mid,:);
            bez.Ce4(reg.sw_mid,:) = newC4(reg.sw_mid,:);     
            %top_slice in elements 1,2,3,4
            bez.Ce1(reg.center_top,:) = newC1(reg.center_top,:);
            bez.Ce1(reg.east_top,:) = newC1(reg.east_top,:);
            bez.Ce1(reg.north_top,:) = newC1(reg.north_top,:);
            bez.Ce1(reg.ne_top,:) = newC1(reg.ne_top,:);        
            bez.Ce2(reg.center_top,:) = newC2(reg.center_top,:);
            bez.Ce2(reg.west_top,:) = newC2(reg.west_top,:);
            bez.Ce2(reg.north_top,:) = newC2(reg.north_top,:);
            bez.Ce2(reg.nw_top,:) = newC2(reg.nw_top,:);
            bez.Ce3(reg.center_top,:) = newC3(reg.center_top,:);
            bez.Ce3(reg.south_top,:) = newC3(reg.south_top,:);
            bez.Ce3(reg.east_top,:) = newC3(reg.east_top,:);
            bez.Ce3(reg.se_top,:) = newC3(reg.se_top,:);                     
            bez.Ce4(reg.center_top,:) = newC4(reg.center_top,:);
            bez.Ce4(reg.south_top,:) = newC4(reg.south_top,:);
            bez.Ce4(reg.west_top,:) = newC4(reg.west_top,:);
            bez.Ce4(reg.sw_top,:) = newC4(reg.sw_top,:);     
            %low_slice in elements 5,6,7,8
            bez.Ce5(reg.center_low,:) = newC5(reg.center_low,:);
            bez.Ce5(reg.east_low,:) = newC5(reg.east_low,:);
            bez.Ce5(reg.north_low,:) = newC5(reg.north_low,:);
            bez.Ce5(reg.ne_low,:) = newC5(reg.ne_low,:);        
            bez.Ce6(reg.center_low,:) = newC6(reg.center_low,:);
            bez.Ce6(reg.west_low,:) = newC6(reg.west_low,:);
            bez.Ce6(reg.north_low,:) = newC6(reg.north_low,:);
            bez.Ce6(reg.nw_low,:) = newC6(reg.nw_low,:);
            bez.Ce7(reg.center_low,:) = newC7(reg.center_low,:);
            bez.Ce7(reg.south_low,:) = newC7(reg.south_low,:);
            bez.Ce7(reg.east_low,:) = newC7(reg.east_low,:);
            bez.Ce7(reg.se_low,:) = newC7(reg.se_low,:);                     
            bez.Ce8(reg.center_low,:) = newC8(reg.center_low,:);
            bez.Ce8(reg.south_low,:) = newC8(reg.south_low,:);
            bez.Ce8(reg.west_low,:) = newC8(reg.west_low,:);
            bez.Ce8(reg.sw_low,:) = newC8(reg.sw_low,:);
            %mid_slice in elements 5,6,7,8
            bez.Ce5(reg.center_mid,:) = newC5(reg.center_mid,:);
            bez.Ce5(reg.east_mid,:) = newC5(reg.east_mid,:);
            bez.Ce5(reg.north_mid,:) = newC5(reg.north_mid,:);
            bez.Ce5(reg.ne_mid,:) = newC5(reg.ne_mid,:);        
            bez.Ce6(reg.center_mid,:) = newC6(reg.center_mid,:);
            bez.Ce6(reg.west_mid,:) = newC6(reg.west_mid,:);
            bez.Ce6(reg.north_mid,:) = newC6(reg.north_mid,:);
            bez.Ce6(reg.nw_mid,:) = newC6(reg.nw_mid,:);
            bez.Ce7(reg.center_mid,:) = newC7(reg.center_mid,:);
            bez.Ce7(reg.south_mid,:) = newC7(reg.south_mid,:);
            bez.Ce7(reg.east_mid,:) = newC7(reg.east_mid,:);
            bez.Ce7(reg.se_mid,:) = newC7(reg.se_mid,:);                     
            bez.Ce8(reg.center_mid,:) = newC8(reg.center_mid,:);
            bez.Ce8(reg.south_mid,:) = newC8(reg.south_mid,:);
            bez.Ce8(reg.west_mid,:) = newC8(reg.west_mid,:);
            bez.Ce8(reg.sw_mid,:) = newC8(reg.sw_mid,:);
        case 8
            bez.Ce5(reg.west_top,:) = newC5(reg.west_top,:);
            bez.Ce5(reg.nw_top,:) = newC5(reg.nw_top,:);
            bez.Ce7(reg.sw_top,:) = newC7(reg.sw_top,:);
            bez.Ce7(reg.west_top,:) = newC7(reg.west_top,:);
        case 9
            bez.Ce1(reg.west_low,:) = newC1(reg.west_low,:);
            bez.Ce1(reg.nw_low,:) = newC1(reg.nw_low,:);
            bez.Ce3(reg.sw_low,:) = newC3(reg.sw_low,:);
            bez.Ce3(reg.west_low,:) = newC3(reg.west_low,:);
        case 10
            bez.Ce6(reg.east_top,:) = newC6(reg.east_top,:);
            bez.Ce6(reg.ne_top,:) = newC6(reg.ne_top,:);
            bez.Ce8(reg.se_top,:) = newC8(reg.se_top,:);
            bez.Ce8(reg.east_top,:) = newC8(reg.east_top,:);
        case 11
            bez.Ce2(reg.east_low,:) = newC2(reg.east_low,:);
            bez.Ce2(reg.ne_low,:) = newC2(reg.ne_low,:);
            bez.Ce4(reg.se_low,:) = newC4(reg.se_low,:);
            bez.Ce4(reg.east_low,:) = newC4(reg.east_low,:);
        case 12
            bez.Ce5(reg.south_top,:) = newC5(reg.south_top,:);
            bez.Ce5(reg.se_top,:) = newC5(reg.se_top,:);
            bez.Ce6(reg.sw_top,:) = newC6(reg.sw_top,:);
            bez.Ce6(reg.south_top,:) = newC6(reg.south_top,:);
        case 13
            bez.Ce1(reg.south_low,:) = newC1(reg.south_low,:);
            bez.Ce1(reg.se_low,:) = newC1(reg.se_low,:);
            bez.Ce2(reg.sw_low,:) = newC2(reg.sw_low,:);
            bez.Ce2(reg.south_low,:) = newC2(reg.south_low,:);
        case 14
            bez.Ce7(reg.north_top,:) = newC7(reg.north_top,:);
            bez.Ce7(reg.ne_top,:) = newC7(reg.ne_top,:);
            bez.Ce8(reg.nw_top,:) = newC8(reg.nw_top,:);
            bez.Ce8(reg.north_top,:) = newC8(reg.north_top,:);
        case 15
            bez.Ce3(reg.north_low,:) = newC3(reg.north_low,:);
            bez.Ce3(reg.ne_low,:) = newC3(reg.ne_low,:);
            bez.Ce4(reg.nw_low,:) = newC4(reg.nw_low,:);
            bez.Ce4(reg.north_low,:) = newC4(reg.north_low,:);
        case 16
            bez.Ce1(reg.sw_mid,:) = newC1(reg.sw_mid,:);
            bez.Ce1(reg.sw_top,:) = newC1(reg.sw_top,:);
            bez.Ce5(reg.sw_low,:) = newC5(reg.sw_low,:);
            bez.Ce5(reg.sw_mid,:) = newC5(reg.sw_mid,:);
        case 17
            bez.Ce2(reg.se_mid,:) = newC2(reg.se_mid,:);
            bez.Ce2(reg.se_top,:) = newC2(reg.se_top,:);
            bez.Ce6(reg.se_low,:) = newC6(reg.se_low,:);
            bez.Ce6(reg.se_mid,:) = newC6(reg.se_mid,:);
        case 18
            bez.Ce3(reg.nw_mid,:) = newC3(reg.nw_mid,:);
            bez.Ce3(reg.nw_top,:) = newC3(reg.nw_top,:);
            bez.Ce7(reg.nw_low,:) = newC7(reg.nw_low,:);
            bez.Ce7(reg.nw_mid,:) = newC7(reg.nw_mid,:);
        case 19
            bez.Ce4(reg.ne_mid,:) = newC4(reg.ne_mid,:);
            bez.Ce4(reg.ne_top,:) = newC4(reg.ne_top,:);
            bez.Ce8(reg.ne_low,:) = newC8(reg.ne_low,:);
            bez.Ce8(reg.ne_mid,:) = newC8(reg.ne_mid,:);                        
    end
end    

%over-write the element_nod indices 
nodes.in1 = newBasisSet{1};
nodes.in2 = newBasisSet{2};
nodes.in3 = newBasisSet{3};
nodes.in4 = newBasisSet{4};
nodes.in5 = newBasisSet{5};
nodes.in6 = newBasisSet{6};
nodes.in7 = newBasisSet{7};
nodes.in8 = newBasisSet{8};

