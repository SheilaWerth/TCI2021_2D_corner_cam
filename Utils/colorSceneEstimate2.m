function [ kHist,v,k ] = colorSceneEstimate2( maxIter,stepsize,penumbraIncident,...
    S,vp,kp,lam2,B,invB,nCoef)
% k is the magnitude of ambient light on the visible side.

m = kp;
u = vp;
ST = S';
kHist = [];
tp = 1;

p = zeros(nCoef,1);
imgDim = sqrt(length(penumbraIncident));

for iter = 1:maxIter
    temp = (m + S*u - penumbraIncident);
    v = u - stepsize*(ST*temp);
    [v,p] = constrained_ell1(v, lam2, B, invB, 100, [0,inf],p);

    k = m - stepsize*(1e4*ones(1,imgDim^2)*temp);
    v(v<0) = 0;
    k(k<0) = 0;
     
     
    t = (sqrt(1+4*tp^2)+1)/2;
    
    u =  v + (tp-1)/t*(v-vp);
    m =  k + (tp-1)/t*(k-kp);

    tp = t;
    vp = v;   
    kp = k;
    
    kHist = [kHist, k];
    cost(iter) = .5*norm(penumbraIncident - k - S*v ,2)^2;

    if mod(iter,100)==1
        figure(77);
        subplot(131)
        plot((1:imgDim)/imgDim*90,v,'.-');
        xlabel('Angle')
        ylabel('Intensity')
        set(gca,'YTick', [])%      
        title('Scene estimate')
        grid on

        subplot(132)
        plot(kHist)
        title('ambient light estimate')
        xlabel('100x Iteration')
        grid on
        drawnow
      
        subplot(133)
        semilogy(cost)
        title('data fidelity cost')
        grid on
        drawnow
        
    end
end

end

