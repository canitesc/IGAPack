function [kk,ff]=feaplyc2sym(kk,ff,bcdof,bcval)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%
%  Synopsis:
%     [kk,ff]=feaplybc(kk,ff,bcdof,bcval)
%
%  Variable Description:
%     kk - system matrix before applying constraints
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value
%
%     For example, there are constraints at d.o.f=2 and 10
%     and their constrained values are 0.0 and 2.5,
%     respectively.  Then, bcdof(1)=2 and bcdof(2)=10; and
%     bcval(1)=1.0 and bcval(2)=2.5.
%-----------------------------------------------------------

bcwt=mean(diag(kk)); % a measure of the average  size of an element in K
ff=ff-kk(:,bcdof)*bcval';  % modify the  force vector
ff(bcdof) = bcwt*bcval;
kk(bcdof,:)=0;  % zero out the rows and  columns of the K matrix
kk(:,bcdof) = 0;
kk(bcdof,bcdof)=bcwt*speye(length(bcdof));  % put ones*bcwt on the diagonal