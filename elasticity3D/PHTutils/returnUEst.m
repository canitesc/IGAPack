function [ up] = returnUEst( p )
%returns the ultra convergent points for the derivative in the reference interval [-1,1]

%define some constants
fudge = 1e-3;
fudge2 = 1e-2;
[~, up3] = quadrature(3, 'GAUSS', 1);
[~, up4] = quadrature(4, 'GAUSS', 1);
[~, up6] = quadrature(6, 'GAUSS', 1);
switch p
    
    case 3
        up = [-1+fudge, 0, 1-fudge];
    case 4
        up = up4';
    case 5
        up = [-1+fudge2, -sqrt(1/3), 0, sqrt(1/3), 1-fudge2];
    case 6
        up = up6';
    case 7
        up = [-1+fudge, -0.769, -0.419, 0, 0.419, 0.769, 1-fudge];
    otherwise
        disp('Invalid p')
        up = [];
end

