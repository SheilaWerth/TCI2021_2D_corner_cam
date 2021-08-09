function [r] = rangeStretch(rVec, indices,N2)
n = length(rVec);
r = ones(N2,1)*NaN;
for i = 1:n
    r(indices(i,1):indices(i,2)) = rVec(i);
end
end

