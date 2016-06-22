function [ u, sigma ] = exactSolHemi( coord )
%exact displacement for the hollow sphere problem
P=1;
b=2;
a=1;
E=1e3;
nu=0.3;

[azimuth,elevation,r] = cart2sph(coord(1),coord(2),coord(3));

phi = azimuth;
theta = pi/2-elevation;


u_r = P*a^3*r/(E*(b^3-a^3))*((1-2*nu)+(1+nu)*b^3/(2*r^3));

[ux, uy, uz] = sph2cart(azimuth,elevation,u_r);
u = [ux, uy, uz];
sigma_r = P*a^3*(b^3-r^3)/(r^3*(a^3-b^3));
sigma_th = P*a^3*(b^3+2*r^3)/(2*r^3*(b^3-a^3));

rot_mat = [sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi);...
    sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi);...
    cos(theta), -sin(theta), 0];

stress_cart = rot_mat*[sigma_r, 0, 0; 0, sigma_th, 0; 0, 0, sigma_th]*rot_mat';

sigma(1) = stress_cart(1,1);
sigma(2) = stress_cart(2,2);
sigma(3) = stress_cart(3,3);
sigma(4) = stress_cart(1,2);
sigma(6) = stress_cart(1,3);
sigma(5) = stress_cart(2,3);


end

