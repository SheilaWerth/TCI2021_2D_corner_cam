function [R,Q,P] = makeRmat_ConstrainedFixJ(sceneRange,knownRange,indices,e,cvec,rvec)
% sceneRange - vector of ranges for each scene element
% knownRange - vector of ranges for each element on floor
% rvec - a vector containing the range of each target in the hidden scene
% cvec - coupling constants
% indices - [startIndexTarget1, stopIndexTarget1; startIndexTarget2, stopIndexTarget2; ...; startIndexTargetN, stopIndexTargetN]

N2 = length(sceneRange);


knownRange = repmat(knownRange,1,N2);
sceneRange = repmat(sceneRange',N2^2,1);

R = (knownRange.^2 + sceneRange.^2 - 2*knownRange.*sceneRange.*e).^-1;
Q = sceneRange.*R;
n = length(cvec);
C = zeros(size(Q));
% loop through targets
for i = 1:n
    R(:,indices(i,1):indices(i,2)) =  rvec(i)*cvec(i)*R(:,indices(i,1):indices(i,2));
    C(:,indices(i,1):indices(i,2)) = cvec(i); 
end
P = C.*(knownRange.^2-sceneRange.^2)./( knownRange.^2 + sceneRange.*(sceneRange-2*knownRange.*e) ).^2;

R(isnan(R)) = 0;

end

