function [ux,uy,uz, stress] = exact_penny(x,E,nu,a,sigma)
% Compute the exact displacement and stresses for the penny-shaped crack using solution
% from Bower book
% Inputs:
%     - x(1,2,3) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu: material properties
%     - sigma : the loading
%     - a: radius of the crack

[theta,r, z]=cart2pol(x(1),x(2), x(3));

rho1 = 1/2*(sqrt((a+r)^2+z^2)-sqrt((a-r)^2+z^2));
rho2 = 1/2*(sqrt((a+r)^2+z^2)+sqrt((a-r)^2+z^2));

drho1dr = (1/4)*(2*a+2*r)/sqrt(a^2+2*a*r+r^2+z^2)-(1/4)*(-2*a+2*r)/sqrt(a^2-2*a*r+r^2+z^2);
drho1dz = (1/2)*z/sqrt(a^2+2*a*r+r^2+z^2)-(1/2)*z/sqrt(a^2-2*a*r+r^2+z^2);
drho2dr = (1/4)*(2*a+2*r)/sqrt(a^2+2*a*r+r^2+z^2)+(1/4)*(-2*a+2*r)/sqrt(a^2-2*a*r+r^2+z^2);
drho2dz = (1/2)*z/sqrt(a^2+2*a*r+r^2+z^2)+(1/2)*z/sqrt(a^2-2*a*r+r^2+z^2);

ur = -nu*sigma*r/E+(1+nu)*sigma*r/(pi*E)*((1-2*nu)*(a*sqrt(rho2^2-a^2)/rho2^2-asin(a/rho2))+2*a^2*abs(z)*sqrt(a^2-rho1^2)/(rho2^2*(rho2^2-rho1^2)));
uz = sigma*z/E +2*sigma*(1+nu)/(pi*E)*(2*(1-nu)*(z/abs(z)*sqrt(a^2-rho1^2)-z*asin(a/rho2))+z*(asin(a/rho2)-a*sqrt(rho2^2-a^2)/(rho2^2-rho1^2)));
%s
ux = ur*cos(theta);
uy = ur*sin(theta);


%compute the derivatives of the displacements with respect to r, z

durdr = -nu*sigma/E+(1+nu)*sigma/pi/E*((1-2*nu)*(a*(rho2^2-a^2)^(1/2)/rho2^2-asin(a/rho2))+2*a^2*abs(z)*(a^2-rho1^2)^(1/2)/rho2^2/(rho2^2-rho1^2))+(1+nu)*sigma*r/pi/E*((1-...
    2*nu)*(a/(rho2^2-a^2)^(1/2)/rho2*drho2dr-2*a*(rho2^2-a^2)^(1/2)/rho2^3*drho2dr+a/rho2^2*drho2dr/(1-a^2/rho2^2)^(1/2))-2*a^2*abs(z)/(a^2-rho1^2)^(1/2)/rho2^2/(rho2^2-...
    rho1^2)*rho1*drho1dr-4*a^2*abs(z)*(a^2-rho1^2)^(1/2)/rho2^3/(rho2^2-rho1^2)*drho2dr-2*a^2*abs(z)*(a^2-rho1^2)^(1/2)/rho2^2/(rho2^2-rho1^2)^2*(2*rho2*drho2dr-2*rho1*drho1dr));

durdz = (1+nu)*sigma*r/pi/E*((1-2*nu)*(a/(rho2^2-a^2)^(1/2)/rho2*drho2dz-2*a*(rho2^2-a^2)^(1/2)/rho2^3*drho2dz+a/rho2^2*drho2dz/(1-a^2/rho2^2)^(1/2))+2*a^2*sign(z)*(a^2-...
    rho1^2)^(1/2)/rho2^2/(rho2^2-rho1^2)-2*a^2*abs(z)/(a^2-rho1^2)^(1/2)/rho2^2/(rho2^2-rho1^2)*rho1*drho1dz-4*a^2*abs(z)*(a^2-rho1^2)^(1/2)/rho2^3/(rho2^2-rho1^2)*drho2dz-...
    2*a^2*abs(z)*(a^2-rho1^2)^(1/2)/rho2^2/(rho2^2-rho1^2)^2*(2*rho2*drho2dz-2*rho1*drho1dz));


duzdz = sigma/E+2*sigma*(1+nu)/pi/E*(2*(1-nu)*(1/abs(z)*(a^2-rho1^2)^(1/2)-z/abs(z)^2*(a^2-rho1^2)^(1/2)*sign(z)-z/abs(z)/(a^2-rho1^2)^(1/2)*rho1*drho1dz-...
    asin(a/rho2)+z*a/rho2^2*drho2dz/(1-a^2/rho2^2)^(1/2))+asin(a/rho2)-a*(rho2^2-a^2)^(1/2)/(rho2^2-rho1^2)+z*(-a/rho2^2*drho2dz/(1-a^2/rho2^2)^(1/2)-...
    a/(rho2^2-a^2)^(1/2)/(rho2^2-rho1^2)*rho2*drho2dz+a*(rho2^2-a^2)^(1/2)/(rho2^2-rho1^2)^2*(2*rho2*drho2dz-2*rho1*drho1dz)));

duzdr = 2*sigma*(1+nu)/pi/E*(2*(1-nu)*(-z/abs(z)/(a^2-rho1^2)^(1/2)*rho1*drho1dr+z*a/rho2^2*drho2dr/(1-a^2/rho2^2)^(1/2))+z*(-a/rho2^2*drho2dr/(1-a^2/rho2^2)^(1/2)-...
    a/(rho2^2-a^2)^(1/2)/(rho2^2-rho1^2)*rho2*drho2dr+a*(rho2^2-a^2)^(1/2)/(rho2^2-rho1^2)^2*(2*rho2*drho2dr-2*rho1*drho1dr)));


epsrr = durdr;
epsthth = ur/r;
epszz = duzdz;
epsrth = 0;
epsthz = 0;
epsrz = 1/2*(durdz+duzdr);

Cmat = E*[(nu-1)/(nu-1+2*nu^2), -nu/(nu-1+2*nu^2), -nu/(nu-1+2*nu^2), 0, 0, 0; -nu/((1+nu)*(-1+2*nu)), (nu-1)/((1+nu)*(-1+2*nu)), -nu/((1+nu)*(-1+2*nu)), 0, 0, 0;...
    -nu/(nu-1+2*nu^2), -nu/(nu-1+2*nu^2), (nu-1)/(nu-1+2*nu^2),0,0,0;0,0,0,1/(2*(1+nu)),0,0;0,0,0,0,1/(2*(1+nu)),0; 0, 0, 0, 0, 0, 1/(2*(1+nu))];

sigma_cyl = Cmat*[epsrr;epsthth; epszz; 2*epsthz; 2*epsrz; 2*epsrth];

trans_mat = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
sigma_cart = trans_mat*[sigma_cyl(1), sigma_cyl(6), sigma_cyl(5); sigma_cyl(6), sigma_cyl(2), sigma_cyl(4); sigma_cyl(5), sigma_cyl(4), sigma_cyl(3)]*trans_mat';
stress = [sigma_cart(1,1); sigma_cart(2,2); sigma_cart(3,3); sigma_cart(1,2); sigma_cart(2,3); sigma_cart(1,3)];


