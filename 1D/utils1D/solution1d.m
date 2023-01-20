function displacements=solution1d(prescribedDof, prescribedValues, stiffness, force, p)
% function to find solution in terms of global displacements

for i=1:length(prescribedDof)
    c=prescribedDof(i);
    stiffness(c,:) = 0;
   % stiffness(:,c) = 0;
    stiffness(c,c) = 1;
    force(c) = prescribedValues(i);    
end
numDof = size(stiffness,2);
disp(['p = ', num2str(p), ' numDof = ', num2str(numDof), ' condition number est. = ', num2str(condest(stiffness))])
displacements = stiffness\force;
