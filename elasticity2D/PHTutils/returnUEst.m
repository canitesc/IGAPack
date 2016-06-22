function [ up] = returnUEst( p )
%returns the ultra convergent points for the derivative in the reference interval [-1,1]

%define some constants
up1 = 1/21*sqrt(147+42*sqrt(7));
up2 = 1/21*sqrt(147-42*sqrt(7));
up3 = 1/33*sqrt(495+66*sqrt(15));
up4 = 1/33*sqrt(495-66*sqrt(15));

switch p
    
    case 3
        up = [-1, 0, 1];
    case 4
        up = [-1, -sqrt(5)/5, sqrt(5)/5, 1];
    case 5
        up = [-1, -sqrt(21)/7, 0, sqrt(21)/7, 1];
    case 6
        up = [-1, -up1, -up2, up2, up1, 1];
    case 7
        up = [-1, -up3, -up4, 0, up4, up3, 1];
    otherwise
        disp('Invalid p')
        up = [];
end

