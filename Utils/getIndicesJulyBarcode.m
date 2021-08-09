function [indices,indicesC,th] = getIndicesJulyBarcode(v,degreeCropStart,degreeCropEnd,thmult)

N2 = length(v);
% indValid1 = round(5/90*N2):round(85/90*N2);
indValid = round(degreeCropStart/90*N2):round(degreeCropEnd/90*N2);
% indValid = round(5/90*N2):round(85/90*N2);

indInvalid = [1:indValid(1), indValid(end):N2];
% th = mean(v(indValid));

indValid2 = round(degreeCropStart/90*N2):round(degreeCropEnd/90*N2);
indInvalid2 = [1:indValid2(1), indValid2(end):N2];

% compute the threshold, from indValid2
th = thmult* mean(v(indValid2));

ind = zeros(N2,1);
% put a one wherever we cross the threshold
ind(v>th) = 1;
% and zero out threshold crossing in invalid region
ind(indInvalid) = 0;
indicesC = 1:N2;
indicesC = indicesC(~logical(ind));

diff_ind = diff(ind);
numObj = sum(diff_ind==1);
indices = zeros(numObj,2);
ind1 = find(diff_ind==1);
ind0 = find(diff_ind==-1);
a = 1:N2;
for ii = 1:numObj
    indices(ii,1) = a(ind1(ii));
    indices(ii,2) = a(ind0(ii));
end
indices = [indices];
