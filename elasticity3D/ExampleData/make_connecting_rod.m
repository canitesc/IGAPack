%script make_connecting_rod.m
%produces the patches for the 3D connecting rod
close all

%big end
addpath ..\nurbs\inst
rad_int = 26.7;
rad_ext = 35.6;
height_conn = 25; %height of the junction with the 2nd patch
center = [0,0,0]; %center of the big end

angle_conn = asin(height_conn/rad_ext);


crv1ext = nrbcirc(rad_ext,center,angle_conn,pi/2);
crv2ext = nrbcirc(rad_ext,center,-angle_conn,angle_conn);
crv3ext = nrbcirc(rad_ext,center,-pi/2,-angle_conn);

crv1int = nrbcirc(rad_int,center,angle_conn,pi/2);
crv2int = nrbcirc(rad_int,center,-angle_conn,angle_conn);
crv3int = nrbcirc(rad_int,center,-pi/2,-angle_conn);


%rotate 90 degrees about the x-axis

rxMat = vecrot(pi/2,[1,0,0]);
crv1ext = nrbtform(crv1ext,rxMat);
crv2ext = nrbtform(crv2ext,rxMat);
crv3ext = nrbtform(crv3ext,rxMat);
crv1int = nrbtform(crv1int,rxMat);
crv2int = nrbtform(crv2int,rxMat);
crv3int = nrbtform(crv3int,rxMat);

%use the Coons method to generate surface NURBS
ptop1 = [0,0,rad_int];
ptop2 = [0,0,rad_ext];
edge1 = nrbline(ptop1,ptop2);
ptopseg1 = [rad_int*cos(angle_conn),0,rad_int*sin(angle_conn)];
ptopseg2 = [rad_ext*cos(angle_conn),0,rad_ext*cos(angle_conn)];
edge2 = nrbline(ptopseg1,ptopseg2);
seg1 = nrbcoons(crv1int, crv1ext, edge1, edge2);

pbotseg1 = [rad_int*cos(angle_conn),0,-rad_int*sin(angle_conn)];
pbotseg2 = [rad_ext*cos(angle_conn),0,-rad_ext*cos(angle_conn)];
edge3 = nrbline(pbotseg1,pbotseg2);
seg2 = nrbcoons(crv2int, crv2ext, edge2, edge3);

pbot1 = [0,0,-rad_int];
pbot2 = [0,0,-rad_ext];
edge4 = nrbline(pbot1, pbot2);
seg3 = nrbcoons(crv3int, crv3ext, edge3, edge4);

%extrude the surfaces to generate the volume NURBS
trans1 = vectrans([0,-11,0]);
trans2 = vectrans([0,-5,0]);
trans3 = vectrans([0,5,0]);

seg1front1 = nrbtform(seg1,trans1);
seg1front2 = nrbtform(seg1,trans2);
seg1front3 = nrbtform(seg1,trans3);

seg2front1 = nrbtform(seg2,trans1);
seg2front2 = nrbtform(seg2,trans2);
seg2front3 = nrbtform(seg2,trans3);

seg3front1 = nrbtform(seg3,trans1);
seg3front2 = nrbtform(seg3,trans2);
seg3front3 = nrbtform(seg3,trans3);

solid1A = nrbextrude(seg1front1, [0,6,0]);
solid1B = nrbextrude(seg1front2, [0,10,0]);
solid1C = nrbextrude(seg1front3, [0,6,0]);

solid2A = nrbextrude(seg2front1, [0,6,0]);
solid2B = nrbextrude(seg2front2, [0,10,0]);
solid2C = nrbextrude(seg2front3, [0,6,0]);

solid3A = nrbextrude(seg3front1, [0,6,0]);
solid3B = nrbextrude(seg3front2, [0,10,0]);
solid3C = nrbextrude(seg3front3, [0,6,0]);


nrbctrlplot(solid1A)
hold on
nrbctrlplot(solid1B)
nrbctrlplot(solid1C)
nrbctrlplot(solid2A)
nrbctrlplot(solid2B)
nrbctrlplot(solid2C)
nrbctrlplot(solid3A)
nrbctrlplot(solid3B)
nrbctrlplot(solid3C)

%input the NURBS object data for the middle patch
coefs4(:,:,1,1) = 1e2*[ -0.050000000000000  -0.050000000000000  -0.050000000000000
   0.253448200000000   0.680001200000000   1.278961000000000
  -0.250000000000000  -0.101235700000000  -0.125000000000000
   0.010000000000000   0.010000000000000   0.010000000000000];
coefs4(:,:,2,1) = [-3.559666000000000  -3.837090500000000  -4.114515000000000
  35.599998966708000  58.251576905151992  98.143855456200001
  -0.000000000218549  -0.000001490132561   0.000000000165784
   0.711933200000000   0.767418100000000   0.822903000000000];
coefs4(:,:,3,1) = 1e2*[-0.050000000000000  -0.050000000000000  -0.050000000000000
   0.253448200000000   0.680001400000000   1.278961000000000
   0.250000000000000   0.101235700000000   0.125000000000000
   0.010000000000000   0.010000000000000   0.010000000000000];
coefs4(:,:,1,2) = 1e2*[0.050000000000000   0.050000000000000   0.050000000000000
   0.253448200000000   0.680001200000000   1.278961000000000
  -0.250000000000000  -0.101235700000000  -0.125000000000000
   0.010000000000000   0.010000000000000   0.010000000000000];
coefs4(:,:,2,2) = [3.559666000000000   3.837090500000000   4.114515000000000
  35.599998966708000  58.251576905151992  98.143855456200001
  -0.000000000218549  -0.000001490132561   0.000000000165784
   0.711933200000000   0.767418100000000   0.822903000000000];
coefs4(:,:,3,2) = 1e2*[0.050000000000000   0.050000000000000   0.050000000000000
   0.253448200000000   0.680001400000000   1.278961000000000
   0.250000000000000   0.101235700000000   0.125000000000000
   0.010000000000000   0.010000000000000   0.010000000000000];

knotU = [0,0,0,1,1,1];
knotV = [0,0,0,1,1,1];
knotW = [0,0,1,1];

solid4 = nrbmak(coefs4([2,1,3,4],:,:,:), {knotU, knotV, knotW})
%solid4.coefs(1,:,:,:) = solid4.coefs(1,:,:,:)+1;
nrbctrlplot(solid4)
hold on

%create the patches for the small end
diam_int = 23;
diam_ext = 44;
height_tail = 25;
center_tail = [146,0,0];


angle_tail = asin(height_tail/diam_ext)
arc1ext = nrbcirc(diam_ext/2,center_tail,pi-angle_tail,pi+angle_tail);
arc2ext = nrbcirc(diam_ext/2,center_tail,pi/2,pi-angle_tail);
arc3ext = nrbcirc(diam_ext/2,center_tail,0,pi/2);
arc4ext = nrbcirc(diam_ext/2,center_tail,-pi/2,0);
arc5ext = nrbcirc(diam_ext/2,center_tail,-pi+angle_tail,-pi/2);

arc1int = nrbcirc(diam_int/2,center_tail,pi-angle_tail,pi+angle_tail);
arc2int = nrbcirc(diam_int/2,center_tail,pi/2,pi-angle_tail);
arc3int = nrbcirc(diam_int/2,center_tail,0,pi/2);
arc4int = nrbcirc(diam_int/2,center_tail,-pi/2,0);
arc5int = nrbcirc(diam_int/2,center_tail,-pi+angle_tail,-pi/2);

%rotate 90 degrees about the x-axis
arc1ext = nrbtform(arc1ext,rxMat);
arc2ext = nrbtform(arc2ext,rxMat);
arc3ext = nrbtform(arc3ext,rxMat);
arc4ext = nrbtform(arc4ext,rxMat);
arc5ext = nrbtform(arc5ext,rxMat);

arc1int = nrbtform(arc1int,rxMat);
arc2int = nrbtform(arc2int,rxMat);
arc3int = nrbtform(arc3int,rxMat);
arc4int = nrbtform(arc4int,rxMat);
arc5int = nrbtform(arc5int,rxMat);


%use the Coons method to generate surface NURBS
qtop1 = center_tail-[diam_int/2*cos(angle_tail),0,-diam_int/2*sin(angle_tail)];
qtop2 = center_tail-[diam_ext/2*cos(angle_tail),0,-diam_ext/2*sin(angle_tail)];
tail_edge2 = nrbline(qtop1,qtop2);
qtopseg1 = center_tail-[diam_int/2*cos(angle_tail),0,diam_int/2*sin(angle_tail)];
qtopseg2 = center_tail-[diam_ext/2*cos(angle_tail),0,diam_ext/2*sin(angle_tail)];
tail_edge1 = nrbline(qtopseg1,qtopseg2);
tail_seg1 = nrbcoons(arc1int, arc1ext, tail_edge1, tail_edge2);

qseg1 = center_tail + [0, 0, diam_int/2];
qseg2 = center_tail + [0, 0, diam_ext/2];
tail_edge3 = nrbline(qseg1,qseg2);
tail_seg2 = nrbcoons(arc2int, arc2ext, tail_edge2, tail_edge3);

qseg3 = center_tail + [diam_int/2, 0, 0];
qseg4 = center_tail + [diam_ext/2, 0, 0];
tail_edge4 = nrbline(qseg3,qseg4);
tail_seg3 = nrbcoons(arc3int, arc3ext, tail_edge3, tail_edge4);

qseg5 = center_tail - [0, 0, diam_int/2];
qseg6 = center_tail - [0, 0, diam_ext/2];
tail_edge5 = nrbline(qseg5,qseg6);
tail_seg4 = nrbcoons(arc4int, arc4ext, tail_edge4, tail_edge5);
tail_seg5 = nrbcoons(arc5int, arc5ext, tail_edge5, tail_edge1);


%extrude the surfaces to generate the volume NURBS

Tseg1front1 = nrbtform(tail_seg1,trans1);
Tseg1front2 = nrbtform(tail_seg1,trans2);
Tseg1front3 = nrbtform(tail_seg1,trans3);

Tseg2front1 = nrbtform(tail_seg2,trans1);
Tseg2front2 = nrbtform(tail_seg2,trans2);
Tseg2front3 = nrbtform(tail_seg2,trans3);

Tseg3front1 = nrbtform(tail_seg3,trans1);
Tseg3front2 = nrbtform(tail_seg3,trans2);
Tseg3front3 = nrbtform(tail_seg3,trans3);

Tseg4front1 = nrbtform(tail_seg4,trans1);
Tseg4front2 = nrbtform(tail_seg4,trans2);
Tseg4front3 = nrbtform(tail_seg4,trans3);

Tseg5front1 = nrbtform(tail_seg5,trans1);
Tseg5front2 = nrbtform(tail_seg5,trans2);
Tseg5front3 = nrbtform(tail_seg5,trans3);


solidT1A = nrbextrude(Tseg1front1, [0,6,0]);
solidT1B = nrbextrude(Tseg1front2, [0,10,0]);
solidT1C = nrbextrude(Tseg1front3, [0,6,0]);

solidT2A = nrbextrude(Tseg2front1, [0,6,0]);
solidT2B = nrbextrude(Tseg2front2, [0,10,0]);
solidT2C = nrbextrude(Tseg2front3, [0,6,0]);

solidT3A = nrbextrude(Tseg3front1, [0,6,0]);
solidT3B = nrbextrude(Tseg3front2, [0,10,0]);
solidT3C = nrbextrude(Tseg3front3, [0,6,0]);

solidT4A = nrbextrude(Tseg4front1, [0,6,0]);
solidT4B = nrbextrude(Tseg4front2, [0,10,0]);
solidT4C = nrbextrude(Tseg4front3, [0,6,0]);

solidT5A = nrbextrude(Tseg5front1, [0,6,0]);
solidT5B = nrbextrude(Tseg5front2, [0,10,0]);
solidT5C = nrbextrude(Tseg5front3, [0,6,0]);


nrbctrlplot(solidT1A)
%hold on
nrbctrlplot(solidT1B)
nrbctrlplot(solidT1C)
nrbctrlplot(solidT2A)
nrbctrlplot(solidT2B)
nrbctrlplot(solidT2C)
nrbctrlplot(solidT3A)
nrbctrlplot(solidT3B)
nrbctrlplot(solidT3C)
nrbctrlplot(solidT4A)
nrbctrlplot(solidT4B)
nrbctrlplot(solidT4C)
nrbctrlplot(solidT5A)
nrbctrlplot(solidT5B)
nrbctrlplot(solidT5C)




% nrbctrlplot(tail_seg1)
% hold on
% nrbctrlplot(tail_seg2)
% nrbctrlplot(tail_seg3)
% nrbctrlplot(tail_seg4)
% nrbctrlplot(tail_seg5)
% 
% 
% 

% nrbctrlplot(arc1ext)
% hold on
% nrbctrlplot(arc2ext)
% nrbctrlplot(arc3ext)
% nrbctrlplot(arc4ext)
% nrbctrlplot(arc5ext)
% 
% nrbctrlplot(arc1int)
% hold on
% nrbctrlplot(arc2int)
% nrbctrlplot(arc3int)
% nrbctrlplot(arc4int)
% nrbctrlplot(arc5int)





save('ConnRod.mat','solid*')



% nrbctrlplot(seg1)
% hold on
% nrbctrlplot(seg2)
% nrbctrlplot(seg3)

% 
% nrbctrlplot(crv1ext)
% hold on
% nrbctrlplot(crv2ext)
% hold on
% nrbctrlplot(crv3ext)
% hold on
% nrbctrlplot(crv1int)
% nrbctrlplot(crv2int)
% nrbctrlplot(crv3int)