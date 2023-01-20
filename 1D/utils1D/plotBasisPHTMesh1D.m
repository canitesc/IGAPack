function  plotBasisPHTMesh1D( PHTelem, GeoMesh, p )

%plots the elements and basis functions stored in PHTelem structure array
%supports multipatches

%we define colors for up to 6 patches
%colorArray = {'blue', 'red', 'green', 'cyan', 'magenta', 'black'};

figure
numPts = 101; %number of plot points to use on each element
uref = linspace(-1,1,numPts);

%1D bernstein polynomials evaluated at the plot points on the master element
B_u = bernstein_basis(uref,p);

numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            %colorIndex = rem((patchIndex-1),6)+1;
            
            umin = PHTelem{patchIndex}(i).vertex(1);
            umax = PHTelem{patchIndex}(i).vertex(2);
            
            coord = zeros(1, numPts);
            R = zeros(numPts, p+1);
            for j=1:numPts
                coord(j) = paramMap1D(GeoMesh{patchIndex}, uref(j), umin, umax);
                tempR = (PHTelem{patchIndex}(i).C)*B_u(j,:)';
                R(j,:) = tempR;
            end            
            
            %plot the p+1 basis functions with support on the current
            %element
            for k=1:p+1
                plot(coord, R(:,k), '-b', 'LineWidth', 2)
                hold on
            end
    %        plot(coord, sum(R')', '-b', 'LineWidth', 2)
    %                    hold on                        

            plot([umin,umax],[0,0],'.r', 'MarkerSize', 10)
            
        end
      
    end
end
axis tight
hold on
drawnow
