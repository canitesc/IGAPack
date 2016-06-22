%script make_c0_cube_with_hole.m
%this requires the nurbs toolbox to be in the path

addpath ../nurbs/inst
close all
clear all
%start with C0 plate with a hole

%initialize geometry on coarsest mesh
rad = 1;
side_fac = 1;

knotU = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
knotV = [0, 0, 1, 1];

coefs(1:3,1,1) = rad.*[-1;0;0];
coefs(1:3,2,1) = rad.*[-0.853553390593274; 0.353553390593274; 0];
coefs(1:3,3,1) = rad.*[-0.603553390593274; 0.603553390593274; 0];
coefs(1:3,4,1) = rad.*[-0.353553390593274;0.853553390593274;0];
coefs(1:3,5,1) = rad.*[0;1;0];

coefs(1:3,1,2) = side_fac.*[-4;0;0];
coefs(1:3,2,2) = side_fac.*[-4;2;0];
coefs(1:3,3,2) = side_fac.*[-4;4;0];
coefs(1:3,4,2) = side_fac.*[-2;4;0];
coefs(1:3,5,2) = side_fac.*[0;4;0];


coefs(4,1,1) = 1;
coefs(4,2,1) = 0.853553390593274;
coefs(4,3,1) = 0.853553390593274;
coefs(4,4,1) = 0.853553390593274;
coefs(4,5,1) = 1;

coefs(4,1,2) = 1;
coefs(4,2,2) = 1;
coefs(4,3,2) = 1;
coefs(4,4,2) = 1;
coefs(4,5,2) = 1;

srf = nrbmak(coefs,{knotU,knotV})
figure
nrbctrlplot(srf)

vol1=nrbrevolve(srf,[0,0,0],[0,1,0],pi/2)
vol1=nrbkntins(vol1,{[0.5,0.5],[],[]})
vol1.coefs(:,3,1,2) = [-4;0;4;1];
% vol1.coefs(:,3,2,2) = [-4;1.707106781186548;4;1];
vol1.coefs(:,3,2,2) = [-3.414213562373095;1.707106781186548;3.414213562373095;0.853553390593274];
vol1.coefs(:,3,3,2) = [-4;4;4;1];
figure
nrbctrlplot(vol1)

%split in 4 patches
coefs = vol1.coefs;
knotU = [0,0,0,1,1,1];
knotV = [0,0,0,1,1,1];
knotW = [0,0,1,1];
coefsSW = coefs(:,1:3,1:3,:);
coefsSE = coefs(:,3:5,1:3,:);
coefsNW = coefs(:,1:3,3:5,:);
coefsNE = coefs(:,3:5,3:5,:);

solidSW = nrbmak(coefsSW,{knotU,knotV,knotW});
solidSE = nrbmak(coefsSE,{knotU,knotV,knotW});
solidNW = nrbmak(coefsNW,{knotU,knotV,knotW});
solidNE = nrbmak(coefsNE,{knotU,knotV,knotW});

figure
nrbctrlplot(solidSW)
hold on
nrbctrlplot(solidSE)
hold on
nrbctrlplot(solidNW)
hold on
nrbctrlplot(solidNE)



%figure
%nrbctrlnetplot3(vol1)
% vol2 = nrbkntins(vol1,{0.5,[],[]})
% % coefs2(:,:,:,1) = vol1.coefs(:,:,:,1);
% % coefs2(:,:,:,3) = vol1.coefs(:,:,:,2);
% % coefs2(:,:,:,2) = (vol1.coefs(:,:,:,1) + vol1.coefs(:,:,:,2))/2;
% % knotW = [0,0,0.5,1,1];
% 
% %vol2 = nrbmak(coefs2,{vol1.knots{1}, vol1.knots{2}, knotW})
% 
% figure
% nrbctrlplot(vol2)
% figure
% nrbplot(vol1,[21,21,21])
% % hold on
% coefs = vol1.coefs;
% nx = size(coefs,2)
% ny = size(coefs,3)
% nz = size(coefs,4)
% for kk=1:nz
%     for jj=1:ny
%         for ii=1:nx
%             plot3(coefs(1,ii,jj,kk)/coefs(4,ii,jj,kk),coefs(2,ii,jj,kk)/coefs(4,ii,jj,kk),coefs(3,ii,jj,kk)/coefs(4,ii,jj,kk), 'r.','MarkerSize',20)
%             [ii, jj, kk]
%             pause
%         end
%     end
% end

            
