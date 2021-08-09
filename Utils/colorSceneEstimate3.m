function [ kHist,v,k ] = colorSceneEstimate3( maxIter,stepsize,penumbraIncident,...
    S,vp,kp,lam2,B,invB,nCoef,lam3,D)
m = kp;
u = vp;
ST = S';
kHist = [];
tp = 1;

p = zeros(nCoef,1);
imgDim = sqrt(length(penumbraIncident));

for iter = 1:maxIter
    temp = (m + S*u - penumbraIncident);
%     v = u - stepsize*(ST*temp + lam2*W_adj(W(u)));
    v = u - stepsize*(ST*temp+lam3*(D'*D)*u);
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
    cost1(iter) = .5*norm(penumbraIncident - k - S*v ,2)^2 ;
    cost2(iter) = lam3*norm(D*v,2)^2;

    if mod(iter,100)==1
        figure(777);
        subplot(141)
        plot((1:imgDim)/imgDim*90,v,'.-');
        set(gca,'YTick', [])
        xlabel('Angle')
        ylabel('Intensity')
        title('Scene estimate')
        grid on

        
        subplot(142)
        plot(kHist)
        title('ambient light estimate')
        xlabel('100x Iteration')
        grid on
        drawnow
        
        
        subplot(143)
        semilogy(cost1)
        title('Data fidelity cost')
        xlabel('100x Iteration')
        grid on
        drawnow
        
        subplot(144)
        semilogy(cost2)
        title('Smoothness cost')
        xlabel('100x Iteration')
        grid on
        drawnow
        
    end
end

end

