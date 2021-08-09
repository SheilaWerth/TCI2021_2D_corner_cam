function [ gC,gR ] = gradRC_fix_knownAmbient( y,k, A,v,f,cVec,...
    indices,R,Q,P)

nT = length(cVec); % number of targets
gC = zeros(nT,1); % gradient of C
gR = zeros(nT,1); % gradient of R

% size(R)
% size(A)
% size(v)
% size(f)
T = ( y - k*f - f.*((R.*A)*v) )';
% t = 1;
% test = ( A(:,indices(t,1):indices(t,2))./Q(:,indices(t,1):indices(t,2)));
% size(T)
% size(test)
% size(( ( A(:,indices(t,1):indices(t,2))./Q(:,indices(t,1):indices(t,2)) )*v ))
for t = 1:nT
    
    gC(t) = T*( -f.*( ( A(:,indices(t,1):indices(t,2)).*Q(:,indices(t,1):indices(t,2)) )*v(indices(t,1):indices(t,2)) ) );
    gR(t) = T*( -f.*( ( A(:,indices(t,1):indices(t,2)).*P(:,indices(t,1):indices(t,2)) )*v(indices(t,1):indices(t,2)) ) );
   
end
end

