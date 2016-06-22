function  plotPHTMesh3DMP( PHTelem, GIFTmesh )
%plots the elements stored in PHTelem structure array
%supports multipatches


%we define colors for up to 6 patches
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};

numPts = 3; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);
wref = linspace(-1,1,numPts);

numPatches = length(PHTelem);


figure
view(-37.50,30)

for patchIndex = 1:numPatches
    
    for i=1:length(PHTelem{patchIndex})
        colorIndex = rem((patchIndex-1),6)+1;
        if isempty(PHTelem{patchIndex}(i).children)
            
            %store the corners of the element in parameter space for easy
            %access
            pxmin = PHTelem{patchIndex}(i).vertex(1);
            pymin = PHTelem{patchIndex}(i).vertex(2);
            pzmin = PHTelem{patchIndex}(i).vertex(3);
            pxmax = PHTelem{patchIndex}(i).vertex(4);
            pymax = PHTelem{patchIndex}(i).vertex(5);
            pzmax = PHTelem{patchIndex}(i).vertex(6);
            
            %determine the edges of the element in the physical space
            coord_up_left = zeros(numPts,3);
            coord_down_left = zeros(numPts,3);
            coord_up_right = zeros(numPts,3);
            coord_down_right = zeros(numPts,3);
            coord_up_front = zeros(numPts,3);
            coord_down_front = zeros(numPts,3);
            coord_up_back = zeros(numPts,3);
            coord_down_back = zeros(numPts,3);
            coord_left_front = zeros(numPts,3);
            coord_right_front = zeros(numPts,3);
            coord_left_back = zeros(numPts,3);
            coord_right_back = zeros(numPts,3);
            
            for j=1:numPts
                coord_up_left(j,:) = paramMap3D(GIFTmesh{patchIndex}, -1, vref(j), 1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_down_left(j,:) = paramMap3D(GIFTmesh{patchIndex}, -1, vref(j), -1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_up_right(j,:) = paramMap3D(GIFTmesh{patchIndex}, 1, vref(j), 1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_down_right(j,:) = paramMap3D(GIFTmesh{patchIndex}, 1, vref(j), -1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_up_front(j,:) = paramMap3D(GIFTmesh{patchIndex}, uref(j), -1, 1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_down_front(j,:) = paramMap3D(GIFTmesh{patchIndex}, uref(j), -1, -1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_up_back(j,:) = paramMap3D(GIFTmesh{patchIndex}, uref(j), 1, 1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_down_back(j,:) = paramMap3D(GIFTmesh{patchIndex}, uref(j), 1, -1, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_left_front(j,:) = paramMap3D(GIFTmesh{patchIndex}, -1, -1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_right_front(j,:) = paramMap3D(GIFTmesh{patchIndex}, 1, -1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_left_back(j,:) = paramMap3D(GIFTmesh{patchIndex}, -1, 1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
                coord_right_back(j,:) = paramMap3D(GIFTmesh{patchIndex}, 1, 1, wref(j), pxmin, pymin, pzmin, pxmax, pymax, pzmax );
            end
            
            %plot the edges of the element
            
            line(coord_up_left(:,1), coord_up_left(:,2), coord_up_left(:,3), 'Color', colorArray{colorIndex})
            line(coord_down_left(:,1), coord_down_left(:,2), coord_down_left(:,3), 'Color', colorArray{colorIndex})
            line(coord_up_right(:,1), coord_up_right(:,2), coord_up_right(:,3), 'Color', colorArray{colorIndex})
            line(coord_down_right(:,1), coord_down_right(:,2), coord_down_right(:,3), 'Color', colorArray{colorIndex})
            line(coord_up_front(:,1), coord_up_front(:,2), coord_up_front(:,3), 'Color', colorArray{colorIndex})
            line(coord_down_front(:,1), coord_down_front(:,2), coord_down_front(:,3), 'Color', colorArray{colorIndex})
            line(coord_up_back(:,1), coord_up_back(:,2), coord_up_back(:,3), 'Color', colorArray{colorIndex})
            
            line(coord_down_back(:,1), coord_down_back(:,2), coord_down_back(:,3), 'Color', colorArray{colorIndex})
            line(coord_left_front(:,1), coord_left_front(:,2), coord_left_front(:,3), 'Color', colorArray{colorIndex})
            line(coord_right_front(:,1), coord_right_front(:,2), coord_right_front(:,3), 'Color', colorArray{colorIndex})
            line(coord_left_back(:,1), coord_left_back(:,2), coord_left_back(:,3), 'Color', colorArray{colorIndex})
            line(coord_right_back(:,1), coord_right_back(:,2), coord_right_back(:,3), 'Color', colorArray{colorIndex})
            
%             %write the element number in the middle
%             coord_mid = paramMap3D(GIFTmesh{patchIndex}, 0, 0, 0, pxmin, pymin, pzmin, pxmax, pymax, pzmax );
%             
%             text(coord_mid(1), coord_mid(2), coord_mid(3), num2str(i), 'Color', colorArray{colorIndex})
        end
    end
end

axis tight
hold on
drawnow
