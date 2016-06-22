function [eldisp,estress]=gbeam(x,y,p,emodu,poisson,d,l)
%-----------------------------------------------------------------------
%  analytical solution:
%  ux=Py/(6EI)*[(6L-3x)*x+(2+poisson)*(y^2-D^2/4)];
%  uy=-P/(6EI)*[3*poisson*y^2*(L-x)+(4+5*poisson)*D^2*x/4+(3*L-x)*x^2];
%-----------------------------------------------------------------------
y=y-d/2;% move (0,0) to below left corner 
inert=d*d*d/12;
pei=p/(6*emodu*inert);
ux=pei*y*((6*l-3*x)*x+(2+poisson)*(y*y-d*d/4));
uy=-pei*(3*poisson*y*y*(l-x)+(4+5*poisson)*d*d*x/4+(3*l-x)*x*x);
dux=pei*6*(l-x)*y;
dvy=-pei*6*poisson*y*(l-x);
dxy=pei*6*(1+poisson)*(y*y-d*d/4);                % shearing strain: duy+dvx

sigx=p*(l-x)*y/inert;
sigy=0;
sigxy=p*(y*y-d*d/4)/2/inert;

eldisp=[ux;uy];
%estrain=[dux;dvy;dxy];
estress=[sigx;sigy;sigxy];

