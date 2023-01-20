function [lhs,rhs] = boundaryCond1D(lhs, rhs, PHTelem, GeoMesh, type)
% function to impose boundary conditions for 1D BVPs
% INPUT: lhs - stiffness (LHS) matrix
%        rhs - RHS vector
%       type - boundary condition types. Possible values 'DD'=Dirichlet at
%       x=0 and x=L; 'DN'=Dirichlet at x=0, Neumann at x=L; 'ND'=Neumann at
%       x=0, Dirichlet at 'x=L'; 'NN'=Neumann at x=0 and x=L
% OUTPUT: lhs - updated lhs matrix
%         rhs - updated rhs matrix
% This function uses exactSol1D to evaluate the exact solution at the
% boundaries

numPatches = length(PHTelem);


switch type
    case 'DD'
        disp('Imposing Dirichlet boundary conditions at x=0 and x=L...')
        bcdof = zeros(1, 2);
        %find the DOFs corresponding to the boundaries
        for patchIndex = 1:numPatches
            for i=1:length(PHTelem{patchIndex})
                if isempty(PHTelem{patchIndex}(i).children)
                    if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                        bcdof(1) = PHTelem{patchIndex}(i).nodes(1);
                    end
                    if (patchIndex==numPatches) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                        bcdof(2) = PHTelem{patchIndex}(i).nodes(end);
                    end
                end
            end
        end
        
        %find the coordinate of the endpoints
        coordLeft = paramMap1D(GeoMesh{1}, -1, 0, 1);
        coordRight = paramMap1D(GeoMesh{numPatches}, 1, 0, 1);
        
        bcval(1) = exact_sol1d(coordLeft);
        bcval(2) = exact_sol1d(coordRight);
        
        [lhs,rhs]=feaplyc2sym(lhs,rhs,bcdof,bcval);
    case 'DN'
        disp('Imposing Dirichlet boundary conditions at x=0 and Neumann at x=L...')
        %find the DOFs corresponding to the boundaries
        for patchIndex = 1:numPatches
            for i=1:length(PHTelem{patchIndex})
                if isempty(PHTelem{patchIndex}(i).children)
                    if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                        bcdof = PHTelem{patchIndex}(i).nodes(1);
                    end
                    if (patchIndex==numPatches) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                        neuNode = PHTelem{patchIndex}(i).nodes(end);
                    end
                end
            end
        end
        
        %find the coordinate of the endpoints
        coordLeft = paramMap1D(GeoMesh{1}, -1, 0, 1);
        bcval = exact_sol1d(coordLeft);
        %bcdof
        
        coordRight = paramMap1D(GeoMesh{numPatches}, 1, 0, 1);
        [~, neuVal] = exact_sol1d(coordRight);
        %neuNode
        rhs(neuNode) = rhs(neuNode) + neuVal;
        %pause
        [lhs,rhs]=feaplyc2sym(lhs,rhs,bcdof,bcval);
    case 'ND'
        disp('Imposing Neumann boundary conditions at x=0 and Dirichlet at x=L...')
        %find the DOFs corresponding to the boundaries
        for patchIndex = 1:numPatches
            for i=1:length(PHTelem{patchIndex})
                if isempty(PHTelem{patchIndex}(i).children)
                    if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                        neuNode = PHTelem{patchIndex}(i).nodes(1);
                    end
                    if (patchIndex==numPatches) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                        bcdof = PHTelem{patchIndex}(i).nodes(end);
                    end
                end
            end
        end
        
        %find the coordinate of the endpoints
        coordLeft = paramMap1D(GeoMesh{1}, -1, 0, 1);
        [~, neuVal] = exact_sol1d(coordLeft);
        rhs(neuNode) = rhs(neuNode) - neuVal;
        
        coordRight = paramMap1D(GeoMesh{numPatches}, 1, 0, 1);
        bcval = exact_sol1d(coordRight);
        
        [lhs,rhs]=feaplyc2sym(lhs,rhs,bcdof,bcval);
        
    case 'NN'
        disp('Imposing Neumann boundary conditions at x=0 and x=L...')
        neuNode = zeros(1, 2);
        %find the DOFs corresponding to the boundaries
        for patchIndex = 1:numPatches
            for i=1:length(PHTelem{patchIndex})
                if isempty(PHTelem{patchIndex}(i).children)
                    if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                        neuNode(1) = PHTelem{patchIndex}(i).nodes(1);
                    end
                    if (patchIndex==numPatches) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                        neuNode(2) = PHTelem{patchIndex}(i).nodes(end);
                    end
                end
            end
        end
        
        %find the coordinate of the endpoints
        coordLeft = paramMap1D(GeoMesh{1}, -1, 0, 1);
        coordRight = paramMap1D(GeoMesh{numPatches}, 1, 0, 1);

        [~, neuVal] = exact_sol1d([coordLeft,coordRight]);        
        rhs(neuNode(1)) = rhs(neuNode(1)) - neuVal(1);
        rhs(neuNode(2)) = rhs(neuNode(2)) + neuVal(2);
                
    otherwise
        error('Case not defined')
end





