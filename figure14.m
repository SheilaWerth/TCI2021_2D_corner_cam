clear; 
clc; clear;close all

warning('off','all')

rmpath(genpath('/Users/sheilawerth/Dropbox/Lab/'))
addpath(genpath('/Users/sheilawerth/Dropbox/Lab/Passive corner camera code review/TCI'))
rmpath(genpath('/Users/sheilawerth/Dropbox/Lab/Passive corner camera code review/TCI/Figure 14/Utils/deleteMe'))
%% PARAMETERS

numRemove = 0; % how many rows should we remove from the top of the penumbra


% load and plot measurement
datStr = 'barcodeReverse_10252019.mat';
load(datStr)

scaleCol = max([test1.r(:)./max(test1.r(:));test1.g(:)./max(test1.g(:));...
    test1.b(:)./max(test1.b(:))]);
colorImage = cat(3,test1.r./max(test1.r(:)),...
    test1.g./max(test1.g(:)),test1.b./max(test1.b(:)))./scaleCol;


figure
imshow(colorImage)
title('Measured Photo')



penumbraIncidentB = test1.b(:)./max(test1.b(:));
penumbraIncidentG = test1.g(:)./max(test1.g(:));
penumbraIncidentR = test1.r(:)./max(test1.r(:));

photoLengthInch = test1.photoLengthInch; % the actual size of the photograph along one edge [inch]


%% LOAD DATA

M = length(penumbraIncidentB); % number of total pixels in the measurement
imgDim = sqrt(M);

pin =  imgDim/photoLengthInch; % pixel per inch, for converting between units

%% A MATRICES AND OTHER TOOLS

load('A_penumbra155.mat') % this is the visibility matrix
A_penumbra = A_penumbra./max(max(A_penumbra));


% camera pixel coordinates
[XX,YY] = meshgrid(1:imgDim,1:imgDim);
XX = fliplr(XX-0.5);
YY = YY-.5;

alpha = atan(XX./YY);
alpha = repmat(alpha(:),[1,imgDim]);

dTheta = (pi/2)/imgDim;
theta = linspace(dTheta/2,pi/2-dTheta/2,imgDim);
theta = repmat(theta,[imgDim^2,1]);

phi = pi/2 + theta + alpha;
e = cos(phi);

% compute range of each camera pixel from the origin (i.e. the occluding
% edge)
knownRange = sqrt(XX.^2 + YY.^2);% assuming each pixel is 1 'unit' long and wide. can convert to inches later...
knownRangeVec = knownRange(:);


% we initially assume the entire scene is in the far field
farFieldRange = 2000; 
% indices is matrix of size N x 2, where N is the number of targets in the
% hidden scene. Entry (n,1) is the starting index (in angle) of the n'th
% target. Entry (n,2) is the last index (in angle) of the n'th target. Our
% model allows us to estimate a single range per each of the N targets.

% we start be assuming everything is in the far field so all angles of the
% hidden scene are at this same, distant, range
indices = [1,155];

% for the R matrix, which describes radial falloff, given scene range
RestScale = makeRmat_ConstrainedFixJ(rangeStretch(farFieldRange,indices,imgDim),knownRangeVec,indices,e,1,farFieldRange);

% A constant used to scale matrices to get magnitudes of unknowns in
% correct order of magnitude
ss = 1e5; 

% put entire forward model together
S =  ss*A_penumbra.*RestScale; 

ST = S';

%% Set up wavelet or dct basis
a = ones(imgDim,1);

level = fix(log2(imgDim));
wavename = 'db4';
[Coef,L] = wavedec(a,level,wavename);
wavebasis = zeros(length(Coef),imgDim);
for ii = 1:imgDim
    x = zeros(imgDim,1);
    x(ii) = 1;
    wavebasis(:,ii) = wavedec(x,level,wavename);
end


B =@(x) wavebasis*x;
invB =@(theta) wavebasis'*theta;



%% Form Initial Estimate of Scene (under far-field assumption)
% v - the 1D hidden scene
% k - ambient light magnitude of visible side
lam2 = 1e-8;


vp = zeros(imgDim,1);
kp = 0;
stepsize = 1e-11;

[ kHist_b,v_b,k_b ] = colorSceneEstimate2( 5000,stepsize,penumbraIncidentB,...
    S,vp,kp,lam2,B,invB,length(Coef));


sceneWhiteFloorB = v_b;
sceneWhiteFloorRefB = v_b;

[ kHist_g,v_g,k_g ] = colorSceneEstimate2( 5000,stepsize,penumbraIncidentG,...
    S,vp,kp,lam2,B,invB,length(Coef));


sceneWhiteFloorG = v_g;
sceneWhiteFloorRefG = v_g;

[ kHist_r,v_r,k_r ] = colorSceneEstimate2( 5000,stepsize,penumbraIncidentR,...
    S,vp,kp,lam2,B,invB,length(Coef));


sceneWhiteFloorR = v_r;
sceneWhiteFloorRefR = v_r;


%% plot initial scene estimates
figure;
subplot(1,2,1)
plot((1:imgDim)/imgDim*90,sceneWhiteFloorB,'.-','Color','b');
hold on
plot((1:imgDim)/imgDim*90,sceneWhiteFloorG,'.-','Color','g');
plot((1:imgDim)/imgDim*90,sceneWhiteFloorR,'.-','Color','r');
grid on
set(gca,'YTick', [])
legend('blue channel','green channel','red channel')
title('Initial scene estimate')

subplot(122)
plot(kHist_b,'Color','b')
hold on
plot(kHist_g,'Color','g')
plot(kHist_r,'Color','r')
title('Ambient light estimates')
legend('blue channel','green channel','red channel')
xlabel('iteration')
grid on
drawnow 

pernumbraWhiteFloor = S*v_b+k_b;

figure;
subplot(1,3,1)
imagesc(reshape(penumbraIncidentB...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('measurement: $y$','Interpreter','latex')
subplot(1,3,2)
imagesc(reshape(pernumbraWhiteFloor...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('$A\hat{v}$ ','Interpreter','latex')
subplot(1,3,3)
imagesc(reshape(pernumbraWhiteFloor - penumbraIncidentB...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('error: $y - A\hat{v}$','Interpreter','latex')
drawnow


%% FIND AND COUNT TARGEST

lowerBoundDeg = 5; % dont perform thresholding at angles less than this
upperBoundDeg = 85; % dont 
thmult =.6 ;%.4

% perform thresholding to count targets and locate targets in angle
[indices1,indicesC,threshold] = getIndicesJulyBarcode((v_b+v_g+v_r)./3,lowerBoundDeg,upperBoundDeg,thmult);
indices = [1 indices1(1,2); (indices1(1,2)+1) (indices1(2,1)-1); indices1(2,1) imgDim ];
numObj = size(indices,1);

% new term to promote smoothness at boundaries between targets
windowWidth = indices(2,2)-indices(2,1)+1;
Dinner = diag(ones(1,windowWidth))+diag(-ones(1,windowWidth-1),1);
Douter = zeros(windowWidth,imgDim);
Douter(:,indices(2,1):indices(2,2)) = Dinner;
Douter = sparse(Douter);


figure;
plot((1:imgDim)/imgDim*90,sceneWhiteFloorB,'.-');
hold on

a = zeros(imgDim,1);
for ii = 1:numObj
    a(indices(ii,1):indices(ii,2)) = threshold;
end
plot((1:imgDim)/imgDim*90,a,'.-');
xlim([0,90])
legend('estimate','threshold')
grid on
xlabel('Angle')
set(gca,'YTick', [])
drawnow 

indicesNotPadded = indices;
[indices] = padIndicesJuly(indices,imgDim);

%% ALTERNATING ESTIMATION

% params
numOuterIter = 40; % number of alternations
stepsize_r = 1;
stepsize_scene = 1e-18;
uB = 1000; % largest allowable range
lB = 90; % smallest allowable range
r = 200*ones(numObj,1); % initialize ranges
rp = r;
u = r;
tp = 1;
maxIterRange = 600;
maxIterScene = 3000;

lam2 = 5e-17;% sparsity
lam3 = 5*1e5*1e8; % smoothness at boundaries



c = 200*ones(3*numObj,1);
d = c;
cp = c;
cHist = [];
rHist = [];
cost = [];

% sceneWhiteFloor = scene/200;

sceneWhiteFloorB = sceneWhiteFloorB/200; % I'm scaling here to have everything that we are estimating around the same order of magnitude
sceneWhiteFloorG = sceneWhiteFloorG/200; 
sceneWhiteFloorR = sceneWhiteFloorR/200; 

[Rest,Qest,Pest] = makeRmat_ConstrainedFixJ(rangeStretch(r,indices,imgDim),knownRangeVec,indices,e,c(1:numObj),r);
Rest_b = Rest; Qest_b = Qest; Pest_b = Pest;
Rest_g = Rest; Qest_g = Qest; Pest_g = Pest;
Rest_r = Rest; Qest_r = Qest; Pest_r = Pest;


%%

for outerIter = 1:numOuterIter
    cHist = [];
    rHist = [];
    cost = [];
    
    
    % RANGE UPDATE
    disp(['Range update ' num2str(outerIter)])
    for iterRange = 1:maxIterRange
        
        
        % compute gradients with respect to unknowns
        [ gC_b,gR_b ] = gradRC_fix_knownAmbient( penumbraIncidentB, k_b, ss*A_penumbra,sceneWhiteFloorB,...
            ones(size(penumbraIncidentB)),c(1:numObj),...
            indices,Rest_b,Qest_b,Pest_b);
        
        [ gC_g,gR_g ] = gradRC_fix_knownAmbient( penumbraIncidentG, k_g, ss*A_penumbra,sceneWhiteFloorG,...
            ones(size(penumbraIncidentG)),c(numObj+1:2*numObj),...
            indices,Rest_g,Qest_g,Pest_g);
        
        [ gC_r,gR_r ] = gradRC_fix_knownAmbient( penumbraIncidentR, k_r, ss*A_penumbra,sceneWhiteFloorR,...
            ones(size(penumbraIncidentR)),c(2*numObj+1:3*numObj),...
            indices,Rest_r,Qest_r,Pest_r);
        
        r = u - stepsize_r*( gR_b + gR_g + gR_r);
        c = d - stepsize_r*( [gC_b; gC_g; gC_r] );
        
        
        r(r<lB) = lB;
        r(r>uB) = uB;

        
        t = (sqrt(1+4*tp^2)+1)/2;
        u =  r + (tp-1)/t*(r-rp);
        d =  c + (tp-1)/t*(c-cp);
        
        tp = t;
        rp = r;
        cp = c;
        
        cHist = [cHist,c];
        rHist = [rHist,r];
        
        
        % update the R matrix with new estimates
        [Rest_b,Qest_b,Pest_b] = makeRmat_ConstrainedFixJ(rangeStretch(r,indices,imgDim),...
            knownRangeVec,indices,e,c(1:numObj),r);
        [Rest_g,Qest_g,Pest_g] = makeRmat_ConstrainedFixJ(rangeStretch(r,indices,imgDim),...
            knownRangeVec,indices,e,c(numObj+1:2*numObj),r);
        [Rest_r,Qest_r,Pest_r] = makeRmat_ConstrainedFixJ(rangeStretch(r,indices,imgDim),...
            knownRangeVec,indices,e,c(2*numObj+1:3*numObj),r);
        
        cost = [cost,1/2*norm(penumbraIncidentB-k_b*ones(size(penumbraIncidentB))...
            -(ss*Rest_b.*A_penumbra)*sceneWhiteFloorB,2)^2 ];
        
        if ~mod(iterRange,100)
        
            figure(100)
            subplot(232)
            
            for ii = 1:numObj
                plot((1/pin)*rHist(ii,:),'Color',(ii/numObj)*[0 1 0])
                hold on; 
            end           
            ylim([(1/pin)*(lB-10),20])
            hold off
            title('Range estimate')
            xlabel('iter')
            ylabel('inch')
            grid on
            legend('obj 1','background','obj 2')
            
                   
            
            subplot(231)
            semilogy(1:length(cost),cost)
            title('data fidelity cost (blue channel)')
            grid on
            xlabel('Iteration')
            
            subplot(234)
            imagesc(reshape(k_b*ones(size(penumbraIncidentB))+...
                (ss*Rest_b.*A_penumbra)*sceneWhiteFloorB,imgDim,imgDim))
            axis image;colorbar;title('$A\hat{v}_{\rm blue}$ ','Interpreter','latex')
            colorbar
            
            subplot(235)
            imagesc(reshape(penumbraIncidentB,imgDim,imgDim))
            axis image;colorbar;title('measurement: $y_{\rm blue}$','Interpreter','latex')
            colorbar
            
            subplot(236)
            imagesc(reshape(penumbraIncidentB - k_b*ones(size(penumbraIncidentB))...
                -(ss*Rest_b.*A_penumbra)*sceneWhiteFloorB,imgDim,imgDim))
            axis image;colorbar;title('error: $y_{\rm blue} - A\hat{v}_{\rm blue}$','Interpreter','latex')
            colorbar
            
            
            
            subplot(233)
            plot(cHist(1,:),'Color','b');
            hold on
            plot(cHist(2*numObj,:),'Color','g')
            plot(cHist(3*numObj,:),'Color','r')           
            grid on
            hold off
            title('Coupling values: c')
            xlabel('iter')
            
          
            sceneScaledB = sceneWhiteFloorB;
            for ii = 1:numObj
                sceneScaledB(indices(ii,1):indices(ii,2)) = ...
                    c(ii) * sceneScaledB(indices(ii,1):indices(ii,2));
            end
            sceneScaledG = sceneWhiteFloorG;
            for ii = 1:numObj
                sceneScaledG(indices(ii,1):indices(ii,2)) = ...
                    c(numObj+ii) * sceneScaledG(indices(ii,1):indices(ii,2));
            end
            
            sceneScaledR = sceneWhiteFloorR;
            for ii = 1:numObj
                sceneScaledR(indices(ii,1):indices(ii,2)) = ...
                    c(2*numObj+ii) * sceneScaledR(indices(ii,1):indices(ii,2));
            end
            

            drawnow
            
        end
    end
    
    
    
    % SCENE UPDATE
    S_b =  ss*A_penumbra.*Rest_b;
    vp_b = sceneWhiteFloorB;
    kp_b = k_b;
    
    S_g =  ss*A_penumbra.*Rest_g;
    vp_g = sceneWhiteFloorG;
    kp_g = k_g;
    
    S_r =  ss*A_penumbra.*Rest_r;
    vp_r = sceneWhiteFloorR;
    kp_r = k_r;
    
    
    disp(['Scene update ' num2str(outerIter)])

    [ kHist_b,v_b,k_b ] = colorSceneEstimate3( maxIterScene,stepsize_scene,penumbraIncidentB,...
        S_b,vp_b,kp_b,lam2,B,invB,length(Coef),lam3,Douter);
    
    [ kHist_g,v_g,k_g ] = colorSceneEstimate3( maxIterScene,stepsize_scene,penumbraIncidentG,...
        S_g,vp_g,kp_g,lam2,B,invB,length(Coef),lam3 ,Douter);
    
    [ kHist_r,v_r,k_r ] = colorSceneEstimate3( maxIterScene,stepsize_scene,penumbraIncidentR,...
        S_r,vp_r,kp_r,lam2,B,invB,length(Coef),lam3 ,Douter);
    
    sceneWhiteFloorB = v_b;
    sceneScaledB = sceneWhiteFloorB;
    
    sceneWhiteFloorG = v_g;
    sceneScaledG = sceneWhiteFloorG;
    
    sceneWhiteFloorR = v_r;
    sceneScaledR = sceneWhiteFloorR;

    


    for ii = 1:numObj
        sceneScaledB(indices(ii,1):indices(ii,2)) = ...
            c(ii) * sceneScaledB(indices(ii,1):indices(ii,2));
    end
    
    for ii = 1:numObj
        sceneScaledG(indices(ii,1):indices(ii,2)) = ...
            c(numObj+ii) * sceneScaledG(indices(ii,1):indices(ii,2));
    end
    
    for ii = 1:numObj
        sceneScaledR(indices(ii,1):indices(ii,2)) = ...
            c(2*numObj+ii) * sceneScaledR(indices(ii,1):indices(ii,2));
    end
    
    plot((1:imgDim)/imgDim*90,sceneScaledB,'.-');
    plot((1:imgDim)/imgDim*90,sceneWhiteFloorRefB,'.-');
    hold off
    grid on
    
    a = zeros(imgDim,1);
    
    
    rp = r;
    u = r;
    
    d = c;
    cp = c;
    tp = 1;
    
    figure(10)
    ntarg = length(r);
    tvec = linspace(0,90,imgDim);
    
    scaleColor = max([sceneScaledR;sceneScaledG; sceneScaledB]);
    for t = 1:ntarg
        l = indices(t,2)-indices(t,1)+1;
        rr = (1/pin)*r(t);
        tt = tvec(indices(t,1):indices(t,2));
        xx = rr*cosd(tt);
        yy = rr*sind(tt);
        
        cplotR = sceneScaledR(indices(t,1):indices(t,2));
        cplotG = sceneScaledG(indices(t,1):indices(t,2));
        cplotB = sceneScaledB(indices(t,1):indices(t,2));
        
        cplotR = cplotR./scaleColor;
        cplotG = cplotG./scaleColor;
        cplotB = cplotB./scaleColor;
        
        scatter(xx,yy,300,[cplotR,cplotG,cplotB],'filled')
        hold on
    end
    % colormap('gray')
    set(gca,'Color','k')
    hold off
    xlim([0 10])
    ylim([0 10])
    axis square
    grid on
    xlabel('X [inch]')
    ylabel('Y [inch]')
    title(['2D reconstruction, Outer iteration = ' num2str(outerIter)])
    drawnow
    
%     
end
 