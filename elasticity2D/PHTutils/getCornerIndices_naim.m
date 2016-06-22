function [ sw, south, se, west, center, east, nw, north, ne ] = getCornerIndices( p,q)
%calculate regions of a (p+1)*(q+1) matrix

alpha = floor((p-1)/2); beta = floor((q-1)/2);

%form the index matrix
%indexMatrix = reshape(1:(p+alpha)*(q+beta),q+alpha,p+beta)';

indexMatrix = reshape(1:(p+1)*(q+1),q+1,p+1)';

%pick the corners
% sw = indexMatrix(1:2,1:2);
% nw = indexMatrix(q:q+1,1:2);
% ne = indexMatrix(q:q+1,p:p+1);
% se = indexMatrix(1:2,p:p+1);
% south = indexMatrix(1:2, 3:p-1);
% west = indexMatrix(3:q-1, 1:2);
% center = indexMatrix(3:q-1, 3:p-1);
% east = indexMatrix(3:q-1,p:p+1);
% north = indexMatrix(q:q+1, 3:p-1);
% % 
% % 
% % %convert to row vectors
% sw = reshape(sw,1,4);
% se = reshape(se,1,4);
% ne = reshape(ne,1,4);
% nw = reshape(nw,1,4);
% north = reshape(north,1,2*(p-3));
% west = reshape(west, 1, 2*(q-3));
% south = reshape(south, 1, 2*(p-3));
% east = reshape(east, 1, 2*(q-3));



%pick the corners
sw = indexMatrix(1:alpha+1,1:alpha+1);
nw = indexMatrix(q-alpha+1:q+1,1:alpha+1);
ne = indexMatrix(q-alpha+1:q+1,p-alpha+1:p+1);
se = indexMatrix(1:alpha+1,p-alpha+1:p+1);
south = indexMatrix(1:alpha+1, (alpha+2):(alpha +2 )*(p-2*alpha-1));
west = indexMatrix((alpha+2):(alpha +2 )*(p-2*alpha-1), 1:alpha+1);
center = indexMatrix((alpha+2):(alpha +2 )*(p-2*alpha-1), (alpha+2):(alpha +2 )*(p-2*alpha-1));
east = indexMatrix((alpha+2):(alpha +2 )*(p-2*alpha-1),q-alpha+1:q+1);
north = indexMatrix(q-alpha+1:q+1, (alpha+2):(alpha +2 )*(p-2*alpha-1));
% 
% 
% %convert to row vectors
sw = reshape(sw,1,(alpha+1)*(beta+1));
se = reshape(se,1,(alpha+1)*(beta+1));
ne = reshape(ne,1,(alpha+1)*(beta+1));
nw = reshape(nw,1,(alpha+1)*(beta+1));
north = reshape(north,1,(alpha+1)*(p-2*alpha-1));
west = reshape(west, 1, (alpha+1)*(p-2*alpha-1));
south = reshape(south, 1, (alpha+1)*(p-2*alpha-1));
east = reshape(east, 1, (alpha+1)*(p-2*alpha-1));


end

