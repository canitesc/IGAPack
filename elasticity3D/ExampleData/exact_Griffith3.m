% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************

function [ux,uy] = exact_Griffith3(x,E,nu,stressState,KI,KII,xTip,adv)
% Compute the exact displacement for the infinite plate with centered crack
% Inputs:
%     - x(1,2) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu, stressState: material properties
%     - sigmato : the loading
%     - xTip,adv,cracklength: crack data

%  inclination of local coord
alfa=atan2(adv(2),adv(1));

% transpose of the rotation matrix
QT=[cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
% local coordinates
xp=QT*(x-xTip)';
% local polar coordinates
[theta,r]=cart2pol(xp(1),xp(2));


%KI=sigmato*sqrt(pi*cracklength);   % exact KI


mu=E/(2*(1+nu));

if ( strcmp(stressState,'PLANE_STRAIN') )
    
    
    kapp = 3-4*nu;
    
    
    %     uxm1 = KI/(2*mu)*sqrt(r/2/pi)*cos(theta/2)*(2-4*nu+2*sin(theta/2)^2);
    %     uym1 = KI/(2*mu)*sqrt(r/2/pi)*sin(theta/2)*(4-4*nu-2*cos(theta/2)^2);
    %
    uxm2 = KII/(2*mu)*sqrt(r/2/pi)*sin(theta/2)*(4-4*nu+2*cos(theta/2)^2);
    %uym2 = KII/(2*mu)*sqrt(r/2/pi)*cos(theta/2)*(-2+4*nu+2*sin(theta/2)^2)
    
    
    uxm1 = KI/(2*mu)*sqrt(r/2/pi)*cos(theta/2)*(kapp - cos(theta));
    uym1 = KI/(2*mu)*sqrt(r/2/pi)*sin(theta/2)*(kapp - cos(theta));
    
    %uxm2 = KII/(2*mu)*sqrt(r/2/pi)*sin(theta/2)*(kapp + 2 + cos(theta));
    uym2 = KII/(2*mu)*sqrt(r/2/pi)*(-1)*cos(theta/2)*(kapp - 2 + cos(theta));
    
    
    %superimpose mode 1 and mode 2 displacements
    ux = uxm1 + uxm2;
    uy = uym1 + uym2;
    
end





