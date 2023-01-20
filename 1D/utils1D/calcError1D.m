function [l2relerr,h1relerr]=calcError1D(PHTelem,GeoMesh,sol0,p)
% calculate the L2 norm and H1 seminorm of the error
numGauss = p+1;
[gaussWeights, gaussLocations] = quadrature(numGauss, 'GAUSS', 1);
gaussLocations=gaussLocations';

l2norm = 0;
h1norm = 0;

l2relerr = 0;
h1relerr = 0;

%define the 1D Bernstein polynomials
[Bu, dBu] = bernstein_basis(gaussLocations, p);
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
            
            for ii=1:numGauss
                
                [ coord, dxdxi]  = paramMap1D( GeoMesh{patchIndex}, gaussLocations(ii), umin, umax );
                [analSol, analDeriv] = exact_sol1d(coord);
                J = det(dxdxi);
                cR = PHTelem{patchIndex}(i).C * Bu(ii,:)';                
                dR = PHTelem{patchIndex}(i).C * dBu(ii,:)'*2/(umax-umin);
                
                % Solve for first derivatives in global coordinates
                dR = dxdxi\dR;                
                
                %calculate displacement values
                compSol = cR'*sol0(sctr);
                
                %calculate the error in stress values
                compDeriv = dR'*sol0(sctr);
                
                l2norm = l2norm + analSol^2*gaussWeights(ii)*J*scalefac;
                h1norm = h1norm + analDeriv^2*gaussWeights(ii)*J*scalefac;
                
                l2relerr = l2relerr + (analSol-compSol)^2*gaussWeights(ii)*J*scalefac;
                h1relerr = h1relerr + (analDeriv-compDeriv)^2*gaussWeights(ii)*J*scalefac;                                                
            end
            
        end
    end
end
l2relerr = sqrt(l2relerr)/sqrt(l2norm);
h1relerr = sqrt(h1relerr)/sqrt(h1norm);