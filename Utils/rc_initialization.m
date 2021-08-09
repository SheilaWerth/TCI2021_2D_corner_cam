function [r,c] = rc_initialization(indices,indices_old,r_old,c_old) 

numObj_old = size(indices_old,1); % previous number of objects
numObj = size(indices,1); % new number of objects

ind_old = cell(numObj_old,1);% initialize a cell array the size of number of previous objects
for ii = 1:numObj_old
    ind_old{ii} = indices_old(ii,1):indices_old(ii,2);
end

% average the old r and c values, because now we have new indices. We will
% have new indices, with no index overlaps to the old ones, initialized as
% the average of all old ones
r = mean(r_old)*ones(numObj,1);
c = mean(c_old)*ones(numObj,1);


% new indices that overlap with old ones will be initialized in that way
for ii = 1:numObj % loop throught the number of current objects
    ind = indices(ii,1):indices(ii,2);
   for jj = 1:numObj_old
       if ~isempty(intersect(ind,ind_old{jj}))
           r(ii) = r_old(jj);
           c(ii) = c_old(jj);
           break
       end
   end
end
