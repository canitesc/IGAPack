function [ B ] = unique_stab( A )
%does the same thing as unique(A, 'stable') with O(n^2) operations
B = [];
for i=1:length(A)
    remove_elem = 0;
    for j=1:i-1
        if A(i)==A(j)
            remove_elem=1;
            break;
        end        
    end
    if remove_elem==0
        B=[B,A(i)];
    end    
end

