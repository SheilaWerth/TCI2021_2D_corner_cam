function [indicesPadded] = padIndicesJuly(indices,N2)
% determine number of objects
nObj = length(indices(:,1));
% iniitalize the new padded indices
indicesPadded = zeros(nObj,2);


ip = 1;

%fill in the padded indices
for p = 1:nObj % loop through the number of objects
    indicesPadded(p,1) = ip;
    if p == nObj
        indicesPadded(p,2) = N2;
    else
        t = floor((indices(p,2)+indices(p+1,1))/2);
        indicesPadded(p,2) = t;
        ip = t+1;
        
    end
end
% indicesPadded(1,2) = 1;
% indicesPadded(2,1) = 2;
end

