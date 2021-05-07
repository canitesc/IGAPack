function [stress] = ghole(x,y,R, tx)
%computes the stress for the plate with a hole problem

[th, r] = cart2pol(x,y);
stressrr = tx/2*(1-R^2/r^2)+tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*th);
stresstt = tx/2*(1+R^2/r^2)-tx/2*(1+3*R^4/r^4)*cos(2*th);
stressrt = -tx/2*(1+2*R^2/r^2-3*R^4/r^4)*sin(2*th);

A = [cos(th)^2, sin(th)^2, 2*sin(th)*cos(th); sin(th)^2, cos(th)^2, -2*sin(th)*cos(th); ...
    -sin(th)*cos(th), sin(th)*cos(th), cos(th)^2-sin(th)^2];

stress = A\[stressrr;stresstt;stressrt];
