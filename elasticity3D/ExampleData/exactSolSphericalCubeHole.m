function [ u, sigma ] = exactSolSphericalCubeHole( coord )
%exact displacement for the hollow sphere problem
S=1;
a=1;
E=1e3;
nu=0.3;
%mu = E/(2*(1+nu));

[azimuth,elevation,r] = cart2sph(coord(1),coord(2),coord(3));

phi = azimuth;
theta = pi/2-elevation;

% A_1 = S*nu/(1+nu);
% A_2 = S*a^5/(7-5*nu);
% A_3 = S*a^3*(6-5*nu)/(2*(7-5*nu));
% B_1 = -S/(2*(1+nu));
% B_2 = -5*S*a^3/(2*(7-5*nu));  
% 
% u_r = 1/(2*mu)*((-A_1*r+(3*A_2)/(2*r^4)-A_3/r^2)+(3*A_1*r-(9*A_2)/(2*r^4)+B_1*(4*nu-2)*r+(B_2*(4*nu-5)/r^2))*cos(theta)*cos(theta));
% 
% u_th = 1/(2*mu)*((-3*A_1*r-(3*A_2)/r^4 + (B_1*r+B_2/r^2)*(2-4*nu))*sin(theta)*cos(theta));
% 
% u_phi = 0;

%disp('first solution')
%[ux, uy, uz] = sph2cart(u_phi,pi/2-u_th,u_r)

%disp('other solution')

ux = (1+nu)*S/(2*E)*((-2*nu)/(1+nu)+(5*nu-6)*a^3/((7-5*nu)*r^3)+(3*a^5)/((7-5*nu)*r^5)*(1-5*coord(3)^2/r^2))*coord(1);
uy = (1+nu)*S/(2*E)*((-2*nu)/(1+nu)+(5*nu-6)*a^3/((7-5*nu)*r^3)+(3*a^5)/((7-5*nu)*r^5)*(1-5*coord(3)^2/r^2))*coord(2);
uz = (1+nu)*S/(2*E)*((2+5*(5-4*nu)*a^3/((7-5*nu)*r^3) + 6*a^5/((7-5*nu)*r^5))*coord(3) +...
    ((-2*nu)/(1+nu)+(5*nu-6)*a^3/((7-5*nu)*r^3)+(3*a^5)/((7-5*nu)*r^5)*(1-5*coord(3)^2/r^2))*coord(3));

u = [ux, uy, uz];
sigma_rr = S*cos(theta)*cos(theta)+S/(7-5*nu)*(a^3/r^3*(6-5*(5-nu)*cos(theta)*cos(theta))+(6*a^5)/r^5*(3*cos(theta)*cos(theta)-1));

sigma_phiphi = 3*S/(2*(7-5*nu))*(a^3/r^3*(5*nu-2+5*(1-2*nu)*cos(theta)*cos(theta))+(a^5)/r^5*(1-5*cos(theta)*cos(theta)));

sigma_thth =  S*sin(theta)*sin(theta)+S/(2*(7-5*nu))*(a^3/r^3*(4-5*nu+5*(1-2*nu)*cos(theta)*cos(theta))+(3*a^5)/r^5*(3-7*cos(theta)*cos(theta)));

sigma_rth =  S*(-1+1/(7-5*nu)*(-5*a^3*(1+nu)/(r^3)+(12*a^5)/r^5))*sin(theta)*cos(theta);


rot_mat = [sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi);...
    sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi);...
    cos(theta), -sin(theta), 0];

stress_cart = rot_mat*[sigma_rr, sigma_rth, 0; sigma_rth, sigma_thth, 0; 0, 0, sigma_phiphi]*rot_mat';

sigma(1) = stress_cart(1,1);
sigma(2) = stress_cart(2,2);
sigma(3) = stress_cart(3,3);
sigma(4) = stress_cart(1,2);
sigma(6) = stress_cart(1,3);
sigma(5) = stress_cart(2,3);


end

