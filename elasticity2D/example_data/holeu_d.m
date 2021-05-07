function [disp] = holeu_d(x, R, E, nu, tx)
%calculates x and y displacement at point x for "Plate with hole" problem

[th, r] = cart2pol(x(1),x(2));

ux = (1+nu)/E*tx*(1/(1+nu)*r*cos(th)+2*R^2/((1+nu)*r)*cos(th)+...
    R^2/(2*r)*cos(3*th)-R^4/(2*r^3)*cos(3*th));
uy = (1+nu)/E*tx*(-nu/(1+nu)*r*sin(th)-(1-nu)*R^2/((1+nu)*r)*sin(th)+...
    R^2/(2*r)*sin(3*th)-R^4/(2*r^3)*sin(3*th));
disp = [ux, uy];