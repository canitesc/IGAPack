function displacements=solution1dNeu(stiffness, force, p)
% function to find solution in terms of global displacements

numDof = size(stiffness,2);
disp(['p = ', num2str(p), ' numDof = ', num2str(numDof), ' condition number est. = ', num2str(condest(stiffness))])
displacements = stiffness\force;
