function [newBasisSet, dimBasis ] = newBasisIndices3D( newBasisVert, elmIn, p, q, r, dimBasis, reg )
%returns the entries of element nod after adding new basis vertices

alpha = floor((p-1)/2);
beta = floor((q-1)/2);
gamma = floor((r-1)/2);

%copy the entries in elmIn to each element in newBasisSet
newBasisSet = cell(1,8);
for i=1:8
    newBasisSet{i} = elmIn;
end

basisCounter = dimBasis;

for j=1:length(newBasisVert)
    switch newBasisVert(j)
        case 1
            %disp('case 1')
            %add the 2*(p-3)*(r-3) deleted basis functions
            newBasisSet{1}(reg.south_mid) = elmIn(reg.south_mid);
            %add the other new basis functions
            newBasisSet{1}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{2}(reg.sw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{2}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{1}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{5}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{1}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{2}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{5}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{6}(reg.sw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{6}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{2}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{6}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
        case 2
            %disp('case 2')
            %add the 2*(q-3)*(r-3) deleted basis functions
            newBasisSet{2}(reg.east_mid) = elmIn(reg.east_mid);
            %add the other new basis functions
            newBasisSet{2}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{4}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{4}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{2}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{6}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{2}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{6}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{8}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{8}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{4}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{8}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 3
            %disp('case 3')
            %add the 2*(p-3)*(r-3) deleted basis functions         
            newBasisSet{3}(reg.north_mid) = elmIn(reg.north_mid);
            %add the other new basis functions
            newBasisSet{3}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{4}(reg.nw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{4}(reg.north_mid) = basisCounter+1:basisCounter+(alpha+1)*(p-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(p-2*beta-1)*(r-2*gamma-1);
            newBasisSet{3}(reg.north_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{7}(reg.north_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter+(alpha+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{3}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.north_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{7}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{8}(reg.nw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{8}(reg.north_mid) = basisCounter+1:basisCounter+(alpha+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{4}(reg.north_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{8}(reg.north_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
        case 4
            %disp('case 4')
            %add the 2*(q-3)*(r-3) deleted basis functions
            newBasisSet{1}(reg.west_mid) = elmIn(reg.west_mid);
            %add the other new basis functions
            newBasisSet{1}(reg.nw_mid) = basisCounter+1:basisCounter+(beta+1)*(alpha+1)*(r-2*gamma-1);
            newBasisSet{3}(reg.sw_mid) = basisCounter+1:basisCounter+(beta+1)*(alpha+1)*(r-2*gamma-1);
            basisCounter = basisCounter+(beta+1)*(alpha+1)*(r-2*gamma-1);
            newBasisSet{3}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{1}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{5}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{1}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{5}(reg.nw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{7}(reg.sw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{7}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{3}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{7}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
        case 5
            %disp('case 5')
            %add 2*(p-3)*(q-3) deleted basis functions
            newBasisSet{1}(reg.center_low) = elmIn(reg.center_low);
            %add the other new basis functions
            newBasisSet{1}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{2}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{2}(reg.center_low) = basisCounter+1:basisCounter+(gamma+1)*(q-2*beta-1)*(p-2*alpha-1);
            basisCounter = basisCounter + (gamma+1)*(q-2*beta-1)*(p-2*alpha-1);
            newBasisSet{1}(reg.north_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{3}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{1}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{2}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.center_low) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{3}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{4}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{4}(reg.center_low) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{2}(reg.north_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{4}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
        case 6
            %disp('case 6')
            %add 2*(p-3)*(q-3) deleted basis functions
            newBasisSet{5}(reg.center_top) = elmIn(reg.center_top);
            %add the other new basis functions
            newBasisSet{5}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{6}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{6}(reg.center_top) = basisCounter+1:basisCounter+(gamma+1)*(q-2*beta-1)*(p-2*alpha-1);
            basisCounter = basisCounter + (gamma+1)*(q-2*beta-1)*(p-2*alpha-1);
            newBasisSet{5}(reg.north_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{7}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{5}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.center_top) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{7}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{8}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{8}(reg.center_top) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter + (gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{6}(reg.north_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{8}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
        case 7
            %disp('case 7')
            %add (p-3)*(q-3)*(r-3) deleted basis functions
            newBasisSet{1}(reg.center_mid) = elmIn(reg.center_mid);
            %add the other new basis functions on the middle slice of
            %elements 1,2,3,4
            newBasisSet{1}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{2}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{2}(reg.center_mid) = basisCounter+1:basisCounter+(q-2*beta-1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (q-2*beta-1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{1}(reg.north_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{3}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{1}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{2}(reg.nw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{3}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{4}(reg.sw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{3}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{3}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(r-2*gamma-1)*(q-2*beta-1);
            newBasisSet{4}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(r-2*gamma-1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(r-2*gamma-1)*(q-2*beta-1);
            newBasisSet{4}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{2}(reg.north_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{4}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);            
            %add the other new basis functions on the top slice of
            %elements 1,2,3,4 and low slice of elements 5,6,7,8
            newBasisSet{1}(reg.center_top) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{5}(reg.center_low) = basisCounter+1:basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            basisCounter = basisCounter+(gamma+1)*(p-2*alpha-1)*(q-2*beta-1);
            newBasisSet{1}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{2}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{5}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{6}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{2}(reg.center_top) = basisCounter+1:basisCounter+(q-2*beta-1)*(p-2*alpha-1)*(beta+1);
            newBasisSet{6}(reg.center_low) = basisCounter+1:basisCounter+(q-2*beta-1)*(p-2*alpha-1)*(beta+1);
            basisCounter = basisCounter + (q-2*beta-1)*(p-2*alpha-1)*(beta+1);
            newBasisSet{1}(reg.north_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{3}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{5}(reg.north_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{7}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{1}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{2}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.center_top) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(alpha+1);
            newBasisSet{7}(reg.center_low) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(alpha+1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(alpha+1);
            newBasisSet{3}(reg.east_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{4}(reg.west_top) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{7}(reg.east_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{8}(reg.west_low) = basisCounter+1:basisCounter+(alpha+1)*(gamma+1)*(q-2*beta-1);
            basisCounter = basisCounter + (alpha+1)*(gamma+1)*(q-2*beta-1);
            newBasisSet{4}(reg.center_top) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(gamma+1);
            newBasisSet{8}(reg.center_low) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(gamma+1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(gamma+1);
            newBasisSet{2}(reg.north_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{4}(reg.south_top) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{6}(reg.north_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            newBasisSet{8}(reg.south_low) = basisCounter+1:basisCounter+(beta+1)*(gamma+1)*(p-2*alpha-1);
            basisCounter = basisCounter + (beta+1)*(gamma+1)*(p-2*alpha-1);
             %add the other new basis functions on the middle slice of
            %elements 5,6,7,8
            newBasisSet{5}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{5}(reg.east_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{6}(reg.west_mid) = basisCounter+1:basisCounter+(alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{6}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{5}(reg.north_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{7}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{5}(reg.ne_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{6}(reg.nw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{7}(reg.se_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{8}(reg.sw_mid) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(r-2*gamma-1);
            basisCounter = basisCounter + (alpha+1)*(beta+1)*(r-2*gamma-1);
            newBasisSet{7}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{7}(reg.east_mid) = basisCounter+1:basisCounter+(r-2*gamma-1)*(q-2*beta-1)*(alpha+1);
            newBasisSet{8}(reg.west_mid) = basisCounter+1:basisCounter+(r-2*gamma-1)*(q-2*beta-1)*(alpha+1);
            basisCounter = basisCounter + (r-2*gamma-1)*(q-2*beta-1)*(alpha+1);
            newBasisSet{8}(reg.center_mid) = basisCounter+1:basisCounter+(p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (p-2*alpha-1)*(q-2*beta-1)*(r-2*gamma-1);
            newBasisSet{6}(reg.north_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            newBasisSet{8}(reg.south_mid) = basisCounter+1:basisCounter+(beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            basisCounter = basisCounter + (beta+1)*(p-2*alpha-1)*(r-2*gamma-1);
            
        case 8
            %add the 4*(q-3) deleted basis functions
            newBasisSet{5}(reg.west_top) = elmIn(reg.west_top);
            %add 8 more basis function near the vertex
            newBasisSet{5}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{7}(reg.west_top) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 9
            %add the 4*(q-3) deleted basis functions
            newBasisSet{1}(reg.west_low) = elmIn(reg.west_low);
            %add 8 more basis function near the vertex
            newBasisSet{1}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{3}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{3}(reg.west_low) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 10
             %add the 4*(q-3) deleted basis functions
            newBasisSet{6}(reg.east_top) = elmIn(reg.east_top);
            %add 8 more basis function near the vertex
            newBasisSet{6}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{8}(reg.east_top) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 11
            %add the 4*(q-3) deleted basis functions
            newBasisSet{2}(reg.east_low) = elmIn(reg.east_low);
            %add 8 more basis function near the vertex
            newBasisSet{2}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{4}(reg.east_low) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 12
            %add the 4*(p-3) deleted basis functions
            newBasisSet{5}(reg.south_top) = elmIn(reg.south_top);
            %add 8 more basis function near the vertex
            newBasisSet{5}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{6}(reg.south_top) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(beta+1)*(gamma+1)*(p-2*alpha-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 13
            %add the 4*(p-3) deleted basis functions
            newBasisSet{1}(reg.south_low) = elmIn(reg.south_low);
            %add 8 more basis function near the vertex
            newBasisSet{1}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{2}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{2}(reg.south_low) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(beta+1)*(gamma+1)*(p-2*alpha-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 14
            %add the 4*(p-3) deleted basis functions
            newBasisSet{7}(reg.north_top) = elmIn(reg.north_top);
            %add 8 more basis function near the vertex
            newBasisSet{7}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{8}(reg.north_top) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(beta+1)*(gamma+1)*(p-2*alpha-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 15
            %add the 4*(p-3) deleted basis functions
            newBasisSet{3}(reg.north_low) = elmIn(reg.north_low);
            %add 8 more basis function near the vertex
            newBasisSet{3}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{4}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{4}(reg.north_low) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(beta+1)*(gamma+1)*(p-2*alpha-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(gamma+1)*(q-2*beta-1);
        case 16
            %add the 4*(r-3) deleted basis functions
            newBasisSet{1}(reg.sw_mid) = elmIn(reg.sw_mid);
            %add 8 more basis function near the vertex
            newBasisSet{1}(reg.sw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{5}(reg.sw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{5}(reg.sw_mid) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1);
        case 17
            %add the 4*(r-3) deleted basis functions
            newBasisSet{2}(reg.se_mid) = elmIn(reg.se_mid);
            %add 8 more basis function near the vertex
            newBasisSet{2}(reg.se_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{6}(reg.se_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{6}(reg.se_mid) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1);
        case 18
            %add the 4*(r-3) deleted basis functions
            newBasisSet{3}(reg.nw_mid) = elmIn(reg.nw_mid);
            %add 8 more basis function near the vertex
            newBasisSet{3}(reg.nw_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{7}(reg.nw_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{7}(reg.nw_mid) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1);
        case 19
            %add the 4*(r-3) deleted basis functions
            newBasisSet{4}(reg.ne_mid) = elmIn(reg.ne_mid);
            %add 8 more basis function near the vertex
            newBasisSet{4}(reg.ne_top) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            newBasisSet{8}(reg.ne_low) = basisCounter+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1);
            %add 4*(q-3) more basis function in the middle of the west_top edge
            newBasisSet{8}(reg.ne_mid) = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+1:basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1); 
            %increment basis counter
            basisCounter = basisCounter+(alpha+1)*(beta+1)*(gamma+1)+(alpha+1)*(beta+1)*(r-2*gamma-1);
    end        
end

%update dimension of basis space
dimBasis = basisCounter;

