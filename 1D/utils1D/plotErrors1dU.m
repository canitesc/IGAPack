function plotErrors1dU( PHTmesh,elementNodes,displacements,C )
%plot the errors the first derivative

numuniqelem = PHTmesh.numberElements;
b_net = PHTmesh.b_net;
p = PHTmesh.p;

numPoints = 51; %number of points to be plotted per element
evalpoints = linspace(-1, 1, numPoints); %evaluation points in the reference element [-1,1]
physPoints = zeros(1, numPoints);
compSol = zeros(1, numPoints);

figure
%plot errors in the first derivative

gaussLocations=linspace(-1,1,2);
numgauss = length(gaussLocations);

compSolg = zeros(numgauss,1);
physPointsg = zeros(numgauss,1);

for k=1:numuniqelem      
    nodes = elementNodes(k,:);
    Ce = C(:,:,k);    
    cpts = b_net(nodes,1);
    w = b_net(nodes,2);
    
    for i=1:numgauss
        pt = gaussLocations(i);
        [~,coordg,~,dNg] = basisFun(pt,cpts,w,Ce,p); 
        compSolg(i) = dNg*displacements(nodes);
        physPointsg(i) = coordg;
    end
    [~,analSolg] = exact_sol1d(physPointsg);
    
    for i=1:numPoints 
        %evaluate the basis functions and calculate physical coordinate and computed solution
        pt = evalpoints(i);
        [~,coord,~,shgradg] = basisFun(pt,cpts,w,Ce,p);        
        compSol(i) = shgradg*displacements(nodes);
        physPoints(i) = coord;
    end
    [~,analSol] = exact_sol1d(physPoints);
       
    plot(physPoints, analSol-compSol, '-r')%, physPointsg, analSolg-compSolg, '.b')
    xaxis = zeros(size(physPoints));
    hold on
    plot( physPoints, xaxis, '-k')
    plot( PHTmesh.knotVector, 0,'.r','MarkerSize',10)
    hold on
end

title('d/dx(u-u_h)')
hold off
drawnow

