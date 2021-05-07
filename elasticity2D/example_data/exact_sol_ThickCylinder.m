function [ disp, stress] =exact_sol_ThickCylinder(x, Emod, nu, rad_int, rad_ext, bound_press)
%exact displacement and stresses for the pressurized cylinder problem

[th, r] = cart2pol(x(1),x(2));


u_r = rad_int^2*bound_press*r/(Emod*(rad_ext^2-rad_int^2))*(1-nu+(rad_ext/r)^2*(1+nu));
disp(1) = u_r*cos(th);
disp(2) = u_r*sin(th);

sigma_rr = rad_int^2*bound_press/(rad_ext^2-rad_int^2)*(1-rad_ext^2/r^2);
sigma_tt = rad_int^2*bound_press/(rad_ext^2-rad_int^2)*(1+rad_ext^2/r^2);
sigma_rt = 0;

A = [cos(th)^2, sin(th)^2, 2*sin(th)*cos(th); sin(th)^2, cos(th)^2, -2*sin(th)*cos(th); ...
    -sin(th)*cos(th), sin(th)*cos(th), cos(th)^2-sin(th)^2];

stress = A\[sigma_rr;sigma_tt;sigma_rt];
