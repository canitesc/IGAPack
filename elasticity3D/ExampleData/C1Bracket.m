%script C1Bracket.m
%generate the geometry for a C1 bracket (multipatch)
clear all
close all

knotU = [0,0,0,1,2,2,2];
knotV = [0,0,1,1];

coefs(1:4, 1, 1) = [0,0,0,1];
coefs(1:4, 2, 1) = [0,1/2,0,1/2];
coefs(1:4, 3, 1) = [1,1/2,0,1/2];
coefs(1:4, 4, 1) = [2,0,0,1];

coefs(1:4, 1, 2) = [-1, 0, 0, 1];
%coefs(1:4, 2, 2) = [-1/2, 1, 0, 1/2];
%coefs(1:4, 3, 2) = [3/2, 1, 0, 1/2];
coefs(1:4, 2, 2) = [-1, 2, 0, 1];
coefs(1:4, 3, 2) = [3, 2, 0, 1];
coefs(1:4, 4, 2) = [3, 0, 0, 1];

% curv = nrbmak(coefs, knotU);
% nrbctrlplot(curv)

surf = nrbmak(coefs, {knotU, knotV});
nrbctrlplot(surf)

figure
%surf2 = nrbdegelev(surf, [0,2]);
%nrbctrlplot(surf2)

knotV1 = [0,0,0,1/2,1/2,1,1,1];
coefsv1(1:4,1,1) = [0, 0, 0, 1];
coefsv1(1:4,2,1) = [-1, 0, 0, 1];
coefsv1(1:4,3,1) = [-1, 0, 0, 1];
coefsv1(1:4,4,1) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefsv1(1:4,5,1) = [-2, 1, 0, 1];

curv1 = nrbmak(coefsv1, knotV1);
%curv1a = nrbdegelev(curv1,1);
nrbctrlplot(curv1)
% figure
% coefsv2 = [ 0       0   0       1
%     -1       0      0       1
%     -1.000000000000000                   0                   0   1.000000000000000
%     -1.000000000000000                   0                   0   1.000000000000000
%   -0.804737854124365   0.471404520791032                   0   0.804737854124365
%   -1.138071187457698   0.804737854124365                   0   0.804737854124365
%   -2.000000000000000   1.000000000000000                   0   1.000000000000000];
% 
% coefsv2 = coefsv2';
% knotV2 = [0, 0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1, 1];
% curv2 = nrbmak(coefsv2, knotV2);
% nrbctrlplot(curv2)
viscircles([-2,0],1)
