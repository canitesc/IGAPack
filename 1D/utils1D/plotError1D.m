function plotError1D( PHTelem, GeoMesh, sol0, p )
%plot the errors the first derivative

% calculate the L2 norm and H1 seminorm of the error
numPlotPts = 101;
plotPts = linspace(-1,1,numPlotPts);


%define the 1D Bernstein polynomials
[Bu, dBu] = bernstein_basis(plotPts, p);
numPatches = length(PHTelem);

for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            umin = PHTelem{patchIndex}(i).vertex(1);
            umax = PHTelem{patchIndex}(i).vertex(2);
            
            %jacobian of the transformation from reference [-1,1]
            %element to the local element in parameter space
            scalefac = (umax - umin)/2;
                        
            nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan
            sctr = PHTelem{patchIndex}(i).nodes(1:nument);
            coord = zeros(1, numPlotPts);
          
            compSol = zeros(1, numPlotPts);
            compDeriv = zeros(1, numPlotPts);
            
            for ii=1:numPlotPts
                
                [ coord(ii), dxdxi]  = paramMap1D( GeoMesh{patchIndex}, plotPts(ii), umin, umax );                                
                
                cR = PHTelem{patchIndex}(i).C * Bu(ii,:)';                
                dR = PHTelem{patchIndex}(i).C * dBu(ii,:)'*2/(umax-umin);
                
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dR;                
                
                %calculate displacement values
                compSol(ii) = cR'*sol0(sctr);
                
                %calculate the error in stress values
                compDeriv(ii) = dR'*sol0(sctr);                               
            end
            
            [analSol, analDeriv] = exact_sol1d(coord);
            
            plot(coord,analDeriv-compDeriv,'-b')            
            hold on
            plot(coord, zeros(1,numPlotPts), '-k')
            hold on
            plot([umin,umax],[0,0],'.r', 'MarkerSize', 10)
            hold on      
            
        end
    end
end
drawnow

