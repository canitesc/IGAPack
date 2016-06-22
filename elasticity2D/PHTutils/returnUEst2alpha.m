function [ up] = returnUEst2alpha( p )
%returns the ultra convergent points for the derivative in the reference interval [-1,1]

%define some constants
%up1 = 1/21*sqrt(147+42*sqrt(7));
%up2 = 1/21*sqrt(147-42*sqrt(7));
%up3 = 1/33*sqrt(495+66*sqrt(15));
%up4 = 1/33*sqrt(495-66*sqrt(15));
fudge = 1e-3;
[~, up4] = quadrature(4, 'GAUSS', 1);
[~, up5] = quadrature(5, 'GAUSS', 1);
[~, up6] = quadrature(6, 'GAUSS', 1);
switch p
    
    case 3
        up = [-1+fudge, 0, 1-fudge];
    case 4
        up = up4';
    case 5
        up = [-1, -sqrt(1/3), 0, sqrt(1/3), 1];
    case 6
        up = [-.792, -0.276, 0.276, 0.792];
    case 7
        up = [-1, -0.528, 0, 0.528, 1];
    otherwise
        disp('Invalid p')
        up = [];
end

