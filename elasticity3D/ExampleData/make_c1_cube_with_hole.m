%script make_c1_cube_with_hole.m
%this requires the nurbs toolbox to be in the path

addpath ../nurbs/inst
close all
clear all
%start with C1 plate with a hole

%initialize geometry on coarsest mesh
rad = 1;
side_fac = 1;

knotU = [0, 0, 0, 0.5, 1, 1, 1];
knotV = [0, 0, 0, 1, 1, 1];

coefs(1:3,1,1) = rad.*[-1;0;0];
coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
coefs(1:3,3,1) = rad.*[-0.353553390593274;0.853553390593274;0];
coefs(1:3,4,1) = rad.*[0;1;0];

coefs(1:3,1,2) = side_fac.*[-2.5,0,0];
coefs(1:3,2,2) = side_fac.*[-2.5,0.75,0];
coefs(1:3,3,2) = side_fac.*[-0.75,2.5,0];
coefs(1:3,4,2) = side_fac.*[0,2.5,0];


coefs(1:3,1,3) = side_fac.*[-4;0;0];
coefs(1:3,2,3) = side_fac.*[-4;4;0];
coefs(1:3,3,3) = side_fac.*[-4;4;0];
coefs(1:3,4,3) = side_fac.*[0;4;0];



coefs(4,1,1) = 1;
coefs(4,2,1) = 0.853553390593274;
coefs(4,3,1) = 0.853553390593274;
coefs(4,4,1) = 1;

coefs(4,1,2) = 1;
coefs(4,2,2) = 1;
coefs(4,3,2) = 1;
coefs(4,4,2) = 1;

coefs(4,1,3) = 1;
coefs(4,2,3) = 1;
coefs(4,3,3) = 1;
coefs(4,4,3) = 1;

coefs

srf = nrbmak(coefs,{knotU,knotV})
figure
nrbctrlplot(srf)
pause
vol1=nrbrevolve(srf,[0,0,0],[0,1,0],pi/2);

%vol1 = nrbpermute (vol1, [2,1,3])
figure
nrbctrlnetplot(vol1)

% 
% coefs2(:,:,:,1) = vol1.coefs(:,:,:,1);
% coefs2(:,:,:,3) = vol1.coefs(:,:,:,2);
% coefs2(:,:,:,2) = (vol1.coefs(:,:,:,1) + vol1.coefs(:,:,:,2))/2;
% knotW = [0,0,0.5,1,1];
% 
% vol2 = nrbmak(coefs2,{vol1.knots{1}, vol1.knots{2}, knotW})
% figure
% nrbctrlplot(vol2)
% 

