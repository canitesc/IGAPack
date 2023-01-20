function [ west, center, east ] = getCornerIndices1D( p )
%calculate regions of a 1*(p+1) matrix
%the regions are
%||----west-----|------center-----|----east-----|| with length:
%||<--alpha+1-->|<--p-2*alpha-1-->|<--alpha+1-->||

alpha = floor((p-1)/2); 

%pick the corners
west = 1:alpha+1;
center = (alpha+2):p-alpha;
east = (p-alpha+1):(p+1);


end

