function [ stressTensor ] = getWedgeStresses( r, theta )



alpha=3*pi/4;
gamma = 2.565819161212361/(2 * alpha);
gamma2=4.281342942588585/(2*alpha);
an=1; bn=0;     % an=1,bn=0 gives pure mode 1 loading and vice versa for mode 2
% %Pass in vectors of r and theta along with the summation term we want
% %(most important is i=1). First order term
% 
% % sigmaxx = r.^( (gamma-1) ) * (gamma) .* (...
% %     an * ( ( 2 + (gamma) + (-1) ) * cos( ((gamma) - 1) * theta) - ( (gamma) - 1 ) .* cos( ((gamma) -3 ) *theta) )...
% %   - bn * ( ( 2 + (gamma) - (-1) ) * sin( ((gamma) - 1) * theta) - ( (gamma) - 1 ) .* sin( ((gamma) -3 ) *theta) ) );
% % 
% % sigmayy = r.^((gamma-1) ) * (gamma) .* (...
% %     an * ( ( 2 - (gamma) - (-1) ) * cos( ((gamma) - 1) * theta) + ( (gamma) - 1 ) .* cos( ((gamma) -3 ) *theta) )...
% %   - bn * ( ( 2 - (gamma) + (-1) ) * sin( ((gamma) - 1) * theta) + ( (gamma) - 1 ) .* sin( ((gamma) -3 ) *theta) ) );
% % 
% % sigmaxy = r.^((gamma-1) ) * (gamma) .* (...
% %     an * ( ( (gamma) - 1 ) * sin( ((gamma) - 3) * theta) - ( (gamma) + (-1) ) .* sin( ((gamma) - 1 ) *theta) )...
% %   + bn * ( ( (gamma) - 1 ) * cos( ((gamma) - 3) * theta) - ( (gamma) - (-1) ) .* cos( ((gamma) - 1 ) *theta) ) );
% 
% 
% 
% sigmarr = an * r.^(gamma - 1) * (gamma + 1) .* ( ...
%           (gamma - 1) * sin((gamma - 1) * alpha) * cos((gamma + 1) * theta) ...
%         - (gamma + 1) * sin((gamma + 1) * alpha) * cos((gamma - 1) * theta) ...
%         - (gamma - 1) * sin((gamma - 1) * alpha) * cos((gamma + 1) * theta) * (gamma + 1) ...
%         + (gamma - 1) * sin((gamma + 1) * alpha) * cos((gamma - 1) * theta) * (gamma - 1) );
%     
% sigmatt = an * r.^(gamma - 1) * gamma * (gamma + 1) .* ( ...
%           (gamma - 1) * sin((gamma - 1) * alpha) * cos((gamma + 1) * theta) ...
%         - (gamma + 1) * sin((gamma + 1) * alpha) * cos((gamma - 1) * theta) );
%     
% sigmart = -an * r.^(gamma - 1) * gamma .* (...
%            -(gamma - 1) * sin((gamma - 1) * alpha) * sin((gamma + 1) * theta) * (gamma + 1)...
%            +(gamma + 1) * sin((gamma + 1) * alpha) * sin((gamma - 1) * theta) * (gamma - 1) );
%        
% sigmaxx = sigmarr .* cos(theta).^2 + sigmatt .* sin(theta).^2 - 2 * sigmart .* sin(theta) .* cos(theta);
% sigmaxy = sigmarr .* sin(theta) .* cos(theta) - sigmatt .* sin(theta) .* cos(theta) + sigmart .* ...
%            (cos(theta).^2 - sin(theta).^2 );
% sigmayy = sigmarr .* sin(theta).^2 + sigmatt .* cos(theta).^2 + 2 * sigmart .* sin(theta) .* cos(theta);


alpha=alpha*2;

Q = -cos((gamma - 1) * alpha/2) / cos((gamma + 1) * alpha/2);

psi_xx1 = (2 - Q * (gamma + 1)) * cos((gamma - 1) * theta) - (gamma - 1) * cos((gamma - 3) * theta);
psi_xy1 = (Q * (gamma + 1)) * sin((gamma - 1) * theta) + (gamma - 1) * sin((gamma - 3) * theta);
psi_yy1 = (2 + Q * (gamma + 1)) * cos((gamma - 1) * theta) + (gamma - 1) * cos((gamma - 3) * theta);

Q2 = -sin((gamma2 - 1) * alpha/2) / sin((gamma2 + 1) * alpha/2);

psi_xx2 = (2 - Q2 * (gamma2 + 1)) * sin((gamma - 1) * theta) - (gamma2 - 1) * sin((gamma2 - 3) * theta);
psi_xy2 = -Q2 * (gamma2 + 1) * cos((gamma2 - 1) * theta) + (gamma2 - 1) * cos((gamma2 - 3) * theta);
psi_yy2 = (2 + Q2 * (gamma2 + 1)) * sin((gamma2 - 1) * theta) + (gamma2 - 1) * sin((gamma2 - 3) * theta);

sigmaxx = an * gamma * r.^(gamma - 1) .* psi_xx1 +  bn * gamma2 * r.^(gamma2 - 1) .* psi_xx2;
sigmaxy = an * gamma * r.^(gamma - 1) .* psi_xy1 +  bn * gamma2 * r.^(gamma2 - 1) .* psi_xy2;
sigmayy = an * gamma * r.^(gamma - 1) .* psi_yy1 +  bn * gamma2 * r.^(gamma2 - 1) .* psi_yy2;

stressTensor = [sigmaxx sigmayy sigmaxy];
% 

end

