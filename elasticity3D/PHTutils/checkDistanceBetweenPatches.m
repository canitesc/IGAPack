%check distance between patches

restoredefaultpath
close all
clear all

addpath ./PHTutils
addpath ./ExampleData
addpath ../nurbs/inst

p=3;
q=3;
r=3;

numBlades = 1;
numPatches = 10*numBlades;
dimBasis = zeros(1, numPatches);
tic
target_rel_error = 1e-3;
targetScale = 0.5;

%dimensions of the domain
L = 1;  %length of the domain
W = 1;  %width of the domain
H = 1;  %height of the domain

dimBasis = zeros(1, numPatches);
%GIFTmesh = init3DGeometryGIFTMP('Propeller',L,W,H,numPatches);
GIFTmesh = init3DGeometryGIFTMP('Propeller',L,W,H,numPatches);
octupleList = cell(numPatches,1);
PHUTelem = cell(numPatches, 1);
numElementsSideX=5;
numElementsSideY=1;
numElementsSideZ=1;

for i=1:numPatches
    
    if round(i/10)==(i/10)
        
        [ PHUTelem{i}, dimBasis(i), quadList{i} ]=initPHTmesh3Dgen( p,q,r, numElementsSideX, numElementsSideY, numElementsSideZ);
        
    else
        
        [PHUTelem{i}, dimBasis(i)] = initPHTmesh3D(p,q,r);
        quadList{i} = 2:9;
        
    end
end

% pause
plotPHTMesh3DMP_1patch( PHUTelem, GIFTmesh,10)

% patchBoundaries = {1,2,4,2;
%     2,3,4,2;
%     1,4,5,6;
%     [2,4],5,[5,4],[6,2];
%     [3,5],6,[5,4],[6,2];
%     4,7,5,6;
%     [5,7],8,[5,4],[6,2];
%     [6,8],9,[5,4],[6,2];
%     5,10,3,4};

%patchBoundaries = {1,2,2,4};

%patchBoundaries ={2,3,2,4};

%patchBoundaries ={1,4,3,1};

%[2,4],5,[3,2],[1,4]

%[3,5],6,[3,2],[1,4]

%patchBoundaries ={4,7,3,1};

%     [5,7],8,[3,2],[1,4];
%     [6,8],9,[3,2],[1,4];

patchBoundaries ={5,10,6,4};
    
patchA= cell2mat(patchBoundaries(1));
patchB=cell2mat(patchBoundaries(2));
faceA=cell2mat(patchBoundaries(3));
faceB=cell2mat(patchBoundaries(4));

plotPHTMesh3DMP_2patches(PHUTelem, GIFTmesh,patchA,patchB) 

plotPHTMesh3DMP_faceColored( PHUTelem, GIFTmesh,patchA,1 )
    
plotPHTMesh3DMP_faceColored( PHUTelem, GIFTmesh,patchB,1 )

figure
plotPHTMesh3DMP_2patches(PHUTelem, GIFTmesh,patchA,patchB)

CheckDistancePhyPlace( PHUTelem, GIFTmesh,patchBoundaries)

checkConformingDist3D_2patches( PHUTelem, GIFTmesh, patchBoundaries, p)