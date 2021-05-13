function [ nrbPatches ] = genPerforatedBeamGeom()
% Generates the geometry for the perforated beam geometry
% Domain description: rectangular beam with lower left corner at (0,-1) and
% upper right corner at (10,1), with circular voids of radius 0.2 centered
% at (2.5, -0.5), (2.5, 0.5), (4,0), (5,-0.5), (5,0.5), (6,0), (7.5, -0.5),
% and (7.5, 0.5)
%  Input: none
%  Output
%  ------
%      nrbPatches : cell array of structures containing the patch geometry
%                   in NURBS toolbox format
%  

nodes = [0, -1; 3.5, -1; 4.5, -1; 5.5, -1; 6.5, -1; 10, -1;
    0, 0; 3.5, 0; 4.5, 0; 5.5, 0; 6.5, 0; 10, 0;
    0, 1; 3.5, 1; 4.5, 1; 5.5, 1; 6.5, 1; 10, 1; 1.5, -1; 1.5, 0; 1.5, 1;
    8.5, -1; 8.5, 0; 8.5, 1];
nodes = [nodes, zeros(size(nodes,1),1)];

% define the quads at the left and right end of the beam
quads = {nodes([1, 19, 20, 7],:), nodes([7, 20, 21, 13],:), ...
           nodes([22, 6, 12, 23],:), nodes([23, 12, 18, 24], :)};


% define the patches corresponding to each of the 8 holes
vertices = {nodes([19,2,8,20],:), nodes([20, 8, 14, 21], :), ...
    nodes([2, 3, 15, 14], :), nodes([3, 4, 10, 9], :), ...
    nodes([9, 10, 16, 15], :), nodes([4, 5, 17, 16], :), ...
    nodes([5, 22, 23, 11], :), nodes([11, 23, 24, 17], :)};

centers = [2.5, -0.5, 0; 2.5, 0.5, 0; 4, 0, 0; 5, -0.5, 0; 5, 0.5, 0;
    6, 0, 0; 7.5, -0.5, 0; 7.5, 0.5, 0];
rad = 0.2;

numPatches = 4*size(centers, 1)+8;
nrbPatches = cell(1, numPatches);
figure
hold on
patchCounter = 0;
% generate the rectangular patches
for i=1:length(quads)
    patchCounter = patchCounter + 1;
    geomData = genQuad(quads{i});
    nrbPatches{patchCounter} = geomData;
end

% generate the patches around the holes
for i=1:size(centers,1)
    geomData = genQuadWInclusion(vertices{i}, centers(i,:), rad);
    if i==3 || i==6
        patchRight = nrbkntins(geomData{2}, {[0.5, 0.5], []});
        patchLeft = nrbkntins(geomData{4}, {[0.5, 0.5], []});
        [nrb_list_right] = splitPatch(patchRight);
        [nrb_list_left] = splitPatch(patchLeft);
        nrbPatches{patchCounter+1} = geomData{1};
        nrbPatches{patchCounter+2} = nrb_list_right{1};
        nrbPatches{patchCounter+3} = nrb_list_right{2};
        nrbPatches{patchCounter+4} = geomData{3};
        nrbPatches{patchCounter+5} = nrb_list_left{1};
        nrbPatches{patchCounter+6} = nrb_list_left{2};
        patchCounter = patchCounter+6;
    else
        for j=1:length(geomData)
            patchCounter = patchCounter+1;
            nrbPatches{patchCounter} = geomData{j};
        end
    end
end

end

