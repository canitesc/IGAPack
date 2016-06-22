function [ stress] = exact_stressesM3( r, theta, KIII  )
% computes the exact stresses in the neighboring of the crack tip
% 
% 
%   stressrr =  1./(4*sqrt(2*pi*r)).*(Kexact1.*(5*cos(theta/2)-cos(3*theta/2))+...
%                   Kexact2.*(-5*sin(theta/2)+3*sin(3*theta/2)))
%   stresstt =  1./(4*sqrt(2*pi*r)).*(Kexact1.*(3*cos(theta/2)+cos(3*theta/2))+...
%                   Kexact2.*(-3*sin(theta/2)-3*sin(3*theta/2)))
%   stressrt =  1./(4*sqrt(2*pi*r)).*(Kexact1.*(sin(theta/2)+sin(3*theta/2))+...
%                   Kexact2.*(cos(theta/2)+cos(3*theta/2)))


% stressrr = 1/sqrt(2*pi*r)*(Kexact1*cos(theta/2)*(1+sin(theta/2)^2)+Kexact2*sin(theta/2)*(1-3*sin(theta/2)^2));
% stresstt = 1/sqrt(2*pi*r)*cos(theta/2)^2*(Kexact1*cos(theta/2)-3*Kexact2*sin(theta/2));
% stressrt = 1/sqrt(2*pi*r)*cos(theta/2)*(Kexact1*sin(theta/2)*cos(theta/2)+Kexact2*(1-3*sin(theta/2)^2));
% 
% 
% 
% A = [cos(theta)^2, sin(theta)^2, 2*sin(theta)*cos(theta); sin(theta)^2, cos(theta)^2, -2*sin(theta)*cos(theta); ...
%     -sin(theta)*cos(theta), sin(theta)*cos(theta), cos(theta)^2-sin(theta)^2];
% stress = A\[stressrr;stresstt;stressrt];
% 
% sxx = stress(1);
% szz = stress(2);
% sxz = stress(3);

sxy = -KIII/sqrt(2*pi*r)*sin(theta/2);
syz = KIII/sqrt(2*pi*r)*cos(theta/2);


stress = [0; 0; 0; sxy; syz; 0];
