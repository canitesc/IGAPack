%script C1Bracket2.m
%generate the geometry for another C1 bracket (multipatch)
clear all
close all

knotU1 = [0,0,0,1/4,1/4,1/2,3/4,3/4,1,1,1];
%knotU1 = [0,0,0,1/6,1/3,1/2,2/3,5/6,1,1,1];
coefsu1(1:4, 1, 1) = [-1,0,0,1];
coefsu2(1:4, 2, 1) = [0,0,0,1];
coefsu1(1:4, 3, 1) = [0,0,0,1];
coefsu1(1:4, 4, 1) = [0,1/2,0,1/2];
coefsu1(1:4, 5, 1) = [1,1/2,0,1/2];
coefsu1(1:4, 6, 1) = [2,0,0,1];
coefsu1(1:4, 7, 1) = [2,0,0,1];
coefsu1(1:4, 8, 1) = [3,0,0,1];


curu1 = nrbmak(coefsu1, knotU1);
nrbctrlplot(curu1)
hold on

knotU2 = [0,0,0,1/4,1/4,1/2,3/4,3/4,1,1,1];

coefsu2(1:4, 1, 1) = [-1,4,0,1];
coefsu2(1:4, 2, 1) = [0,4,0,1];
coefsu2(1:4, 3, 1) = [0,4,0,1];
coefsu2(1:4, 4, 1) = [0,3/2,0,1/2];
coefsu2(1:4, 5, 1) = [1,3/2,0,1/2];
coefsu2(1:4, 6, 1) = [2,4,0,1];
coefsu2(1:4, 7, 1) = [2,4,0,1];
coefsu2(1:4, 8, 1) = [3,4,0,1];



curu2 = nrbmak(coefsu2, knotU2);
nrbctrlplot(curu2)
hold on


knotV1 = [0,0,1,1];
coefsv1(1:4,1,1) = [-1, 0, 0, 1];
coefsv1(1:4,2,1) = [-1, 4, 0, 1]; 

curv1 = nrbmak(coefsv1, knotV1);
%curv1a = nrbdegelev(curv1,1);
nrbctrlplot(curv1)


knotV2 = [0,0,1,1];
coefsv2(1:4,1,1) = [3, 0, 0, 1];
coefsv2(1:4,2,1) = [3, 4, 0, 1]; 

curv2 = nrbmak(coefsv2, knotV2);
%curv1a = nrbdegelev(curv1,1);
hold on
nrbctrlplot(curv2)

srf = nrbcoons(curu1, curu2, curv1, curv2);
figure
nrbctrlplot(srf)
