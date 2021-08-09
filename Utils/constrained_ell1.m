function [fhat,p,diffP] = constrained_ell1(y, lam, W, invW, maxIter, bounds, p)

lowerB = bounds(1);
upperB = bounds(2);

diffP = zeros(maxIter,1);
q = p;
tp = 1;
stepsize = 1/lam^2;
for iter = 1:maxIter
    r = y - lam*invW(q);
    r(r<lowerB) = lowerB;
    r(r>upperB) = upperB;
    grad = -lam*W(r);

    pnext = q - stepsize*grad;
    
    % the shrinkage operator
    pnext(pnext<-1) = -1;
    pnext(pnext>1) = 1;
    
    t = (sqrt(1+4*tp^2)+1)/2;
    q = pnext + (tp-1)/t*(pnext - p);
    
    diffP(iter) = norm(pnext-p);
%     if diffP(iter) < 1e-10
%         break
%     end
    tp = t;
    p = pnext;
end

fhat = y - lam*invW(p);
fhat(fhat<lowerB) = lowerB;
fhat(fhat>upperB) = upperB;