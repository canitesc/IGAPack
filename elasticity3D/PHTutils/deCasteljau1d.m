function [ row1, row2 ] = deCasteljau1d( row )
% Split row array into two arrays of the same size using deCasteljau
% algorithm ( http://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm )

num_ent = size(row,2);
row2 = zeros(1, num_ent);
row1 = zeros(1, num_ent);
temp = zeros(1, num_ent);

temp(1,:) = row;
row1(1) = row(1);
row2(end) = row(end);
for j=1:num_ent-1
    for k=1:num_ent-j
        temp(k) = (temp(k)+temp(k+1))/2;
    end    
    row1(j+1) = temp(1);
    row2(end-j) = temp(num_ent-j);
end
