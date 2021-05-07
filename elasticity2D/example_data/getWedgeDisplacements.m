function [ u, v ] = getWedgeDisplacements( r, theta, shearMod, kappa )
% Pass in vectors of r and theta along with the summation term we want
% (most important is i=1). First order term

% u = r.^((gamma)) / (2*shearMod) .* (...
%     an * ( ( kappa + (gamma) + (-1) ) * cos((gamma) * theta) - (gamma) .* cos( ((gamma) -2 ) *theta) )...
%   - bn * ( ( kappa + (gamma2) - (-1) ) * sin((gamma2) * theta) - (gamma2) .* sin( ((gamma2) -2 ) *theta) ) );
% 
% v = r.^((gamma)) / (2*shearMod) .* (...
%     an * ( ( kappa - (gamma) - (-1) ) * sin((gamma) * theta) + (gamma) .* sin( ((gamma) -2 ) *theta) )...
%   + bn * ( ( kappa - (gamma2) + (-1) ) * cos((gamma2) * theta) + (gamma2) .* cos( ((gamma2) -2 ) *theta) ) );

% calculate the displacements for the wedge problem. These are taken
% from the dominant terms of the Williams exansion. 
% 

alpha=3*pi/4;
gamma = 2.565819161212361/(2 * alpha);
gamma2=4.281342942588585/(2*alpha);
an=1; bn=0;     % an=1,bn=0 gives pure mode 1 loading and vice versa for mode 2


alpha = alpha * 2;
Q = -cos((gamma - 1) * alpha/2) / cos((gamma + 1) * alpha/2);

psi_1x = 1/(2*shearMod) * ( ( kappa - Q * (gamma + 1) )  * cos(gamma * theta) ...
    - gamma * cos((gamma - 2) * theta));
psi_1y = 1/(2*shearMod) * ( ( kappa + Q * (gamma + 1) )  * sin(gamma * theta) ...
    + gamma * sin((gamma - 2) * theta));

Q2 = -sin((gamma2 - 1) * alpha/2) / sin((gamma2 + 1) * alpha/2);

psi_2x = 1/(2*shearMod) * ( ( kappa - Q2 * (gamma2 + 1) )  * sin(gamma2 * theta) ...
    - gamma2 * sin((gamma2 - 2) * theta));
psi_2y = 1/(2*shearMod) * ( -( kappa + Q2 * (gamma2 + 1) )  * cos(gamma2 * theta) ...
    - gamma2 * cos((gamma2 - 2) * theta));

u = an * r.^gamma .* psi_1x + bn * r.^gamma .* psi_2x;
v = an * r.^gamma .* psi_1y + bn * r.^gamma .* psi_2y;

end

