function [nurbs_list] = splitPatch(nurbs)
% splits a patch that contains C^0 lines. This performs one splitting, for
% multiple C^0 lines, it should be applied recursively
% Input:
% ------
%     nurbs - structure containing the geometry information in the NURBS
%             toolbox format.
% Output:
% -------
%     nurbs_list : list of structures which cover the same geometry but 
%                  without C^0 lines

p = nurbs.order(1)-1;
q = nurbs.order(2)-1;

% get the knot multiplicities
[uniqueKnotU, multU] = getMultiplicity(nurbs.knots{1});
[uniqueKnotV, multV] = getMultiplicity(nurbs.knots{2});

% trim the first and the last entries in the knot vectors
trimmedKnotU = uniqueKnotU(2:end-1);
multU = multU(2:end-1);
trimmedKnotV = uniqueKnotV(2:end-1);
multV = multV(2:end-1);

% loop over the knots in the u direction to find the index of the C^0 line
for i=1:length(trimmedKnotU)
    if multU(i)==p
        disp(['Found C^0 line at u=', num2str(trimmedKnotU(i))])
        [~,j] = find(nurbs.knots{1}==trimmedKnotU(i));
        ic0 = j(1)-1;
        ic1 = j(end)-1;
        coefs1 = nurbs.coefs(:,1:ic0,:);
        coefs2 = nurbs.coefs(:,ic0:end,:);
        knotU1 = normalizeKnot(nurbs.knots{1}([1:ic1+1,ic1+1]));
        knotU2 = normalizeKnot(nurbs.knots{1}([ic0+1,ic0+1:end]));
        patch1 = nrbmak(coefs1, {knotU1, nurbs.knots{2}});
        patch2 = nrbmak(coefs2, {knotU2, nurbs.knots{2}});
        nurbs_list = {patch1, patch2};
        return
    end   
end

% loop over the knots in the v direction to find the index of the C^0 line
for i=1:length(trimmedKnotV)
    if multV(i)==q
        disp(['Found C^0 line at v=', num2str(trimmedKnotV(i))])
        [~,j] = find(nurbs.knots{2}==trimmedKnotV(i));
        ic0 = j(1);
        coefs1 = nurbs.coefs(:,:,1:ic0);
        coefs2 = nurbs.coefs(:,:,ic0:end);
        knotV1 = normalizeKnot(nurbs.knots{2}([1:ic0,ic0]));
        knotV2 = normalizeKnot(nurbs.knots{2}([ic0,ic0:end]));
        patch1 = nrbmak(coefs1, {nurbs.knots{1}, knotV1});
        patch2 = nrbmak(coefs2, {nurbs.knots{2}, knotV2});
        nurbs_list = {patch1, patch2};
        return
    end   
end
