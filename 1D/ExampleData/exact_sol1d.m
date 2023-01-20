function [ y, dy, f, ddy ] = exact_sol1d( x )
% compute the exact solution and derivatives

[ a0, a1 ] = coefPDE( x );
% y = exp(x);
% dy = exp(x);
% ddy = exp(x);
% f = -exp(x);
% 
% 
% y = x.^2;
% dy = 2.*x;
% ddy = 2;
% f = -ddy;
% % 
% 
% y = x-1;
% dy = ones(size(x));
% ddy = 0;
% f = -ddy;

% 


%  y = x.^3;
%  dy = 3.*x.^2;
%  ddy = 6.*x;
% f = -ddy;

%   y = x.^4;
%   dy = 4.*x.^3;
%   ddy = 12.*x.^2;
%  f = -ddy;
% 
% y = x.^5;
% dy = 5.*x.^4;
% ddy = 20.*x.^3;
% f = -ddy;

% y = x.^6;
% dy = 6.*x.^5;
% ddy = 30.*x.^4;
% f = -ddy;

% y = x.^7;
% dy = 7.*x.^6;
% ddy = 42.*x.^5;
% f = -ddy;
% 
% y = x.^8;
% dy = 8.*x.^7;
% ddy = 56.*x.^6;
% f = -ddy;

% 
% 

Pe=a1;
c1 = 1/(1+exp(Pe));
y = -c1+c1*exp(Pe*x);
dy = c1*Pe*exp(Pe*x);
ddy = c1*Pe^2*exp(Pe*x);

% 
% y = x.^(3/2);
% dy = 3/2.*x.^(1/2);
% ddy = 3/4.*x.^(-1/2);
% k=a0;
% y = sin(k*x);
% dy = k*cos(k*x);
% ddy = -k^2*sin(k*x);
f = -ddy + a1*dy + a0*y;
