function [ PHTelem, scaleCoef ] = scaleBasis(PHTelem, dimBasis)
%Scales the Bezier coefficients for the basis defined by PHTelem to improve
%conditioning

scaleCoef = zeros(1, dimBasis);


%loop over the active elements and find the maximum Bezier ordinate
%corresponding to each basis function
for e=1:length(PHTelem)
    if isempty(PHTelem(e).children)
        nument = size(PHTelem(e).C, 1);
        nodes = PHTelem(e).nodes(1:nument);        
        maxVals = max(PHTelem(e).C,[],2);                
        scaleCoef(nodes) = max([scaleCoef(nodes);maxVals']);        
    end
end

%we loop over the active elements and scale the basis functions by the
%reciprocal of the maximum
scaleCoef = 1./scaleCoef;
for e=1:length(PHTelem)
    if isempty(PHTelem(e).children)
        nument = size(PHTelem(e).C, 1);
        nodes = PHTelem(e).nodes(1:nument);
        for i=1:nument
            PHTelem(e).C(i,:) = PHTelem(e).C(i,:)*scaleCoef(nodes(i));
        end
    end
end
        


end

