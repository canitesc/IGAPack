function  plotPHTMeshMP( PHTelem, GIFTmesh )
%plots the elements stored in PHTelem structure array
%supports multipatches

%we define colors for up to 6 patches
colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};

numPts = 11; %number of plot points to use on each edge
uref = linspace(-1,1,numPts);
vref = linspace(-1,1,numPts);

numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            colorIndex = rem((patchIndex-1),6)+1;
            %store the corners of the element in parameter space for easy
            %access
            pxmin = PHTelem{patchIndex}(i).vertex(1);
            pymin = PHTelem{patchIndex}(i).vertex(2);
            pxmax = PHTelem{patchIndex}(i).vertex(3);
            pymax = PHTelem{patchIndex}(i).vertex(4);
            
            %determine the edges of the element in the physical space
            coord_south = zeros(numPts,2);
            coord_east = zeros(numPts,2);
            coord_north = zeros(numPts,2);
            coord_west = zeros(numPts,2);
            
            
            for j=1:numPts
                coord_south(j,:) = paramMap(GIFTmesh{patchIndex}, uref(j), -1, pxmin, pymin, pxmax, pymax );
                coord_east(j,:) = paramMap(GIFTmesh{patchIndex}, 1, vref(j), pxmin, pymin, pxmax, pymax );
                coord_north(j,:) = paramMap(GIFTmesh{patchIndex}, -uref(j), 1, pxmin, pymin, pxmax, pymax );
                coord_west(j,:) = paramMap(GIFTmesh{patchIndex}, -1, -vref(j), pxmin, pymin, pxmax, pymax );
            end
            
            %plot the edges of the element
            coords = [coord_south; coord_east; coord_north; coord_west];
            line(coords(:,1), coords(:,2), 'Color', colorArray{colorIndex})
            
%             %write the element number in the middle
%             coord_mid = paramMap(GIFTmesh{patchIndex}, 0, 0, pxmin, pymin, pxmax, pymax );
%             
%             text(coord_mid(1), coord_mid(2), num2str(i), 'Color', colorArray{patchIndex})
        end
    end
end
axis tight
hold on
drawnow
