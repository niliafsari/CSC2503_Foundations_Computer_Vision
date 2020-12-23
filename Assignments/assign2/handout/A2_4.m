%% File: dinoTestF
%% A2 2017 handout code
%% Estimate F matrix from synthetic corresponding points.
%%
%% ADJ

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;

global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  cd ~jepson/pub/matlab   
  startup;
  cd(dir);
end

reconRoot = '.';  %% CHANGE THIS to the directory you installed A4 handout/
addpath([reconRoot '/utils']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  %% Number of ransac trials to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up cameras
%% The cameras are automatically rotated by projectDino to fixate
%% on the mean of the 3D points.  We do not always want to allow
%% the cameras to move in this fashion (eg. when we change sclZ).
%% So we will compute the rotations for the left and right cameras
%% once and for all, and then use these.
f = 100; % focal length
dLeft = [-50, 0, -150]';  % Center of projection for left camera
dRight = [50, 0, -150]';  % Center of projection for right camera
%% Compute camera rotations to fixate on Dino's center.


[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 1);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 1);
Rright = MextRight(:, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data...
sclZ = 1;
%% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

%% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);

U=[MextRight(1:3,1:3)' -MextRight(1:3,1:3)'*MextRight(1:3,4) ; 0 0 0 1 ];
Rt=MextLeft*U;
R=Rt(1:3,1:3);
tau=Rt(1:3,4);
tau_t=[0 -tau(3) tau(2);tau(3) 0 -tau(1); -tau(2)  tau(1) 0];
E=tau_t*R;
F_0=inv(MintLeft)*E*inv(MintRight);

% Build correspondence data
clear imPts;
Ntry=25;
sigma1=logspace(-5,1.65,80);
err_sig=zeros(1,length(sigma1));
err_sig1=zeros(1,length(sigma1));
for v=1:length(sigma1)
    v
    maxperpErr=zeros(1,Ntry);
    maxperpErr1=zeros(1,Ntry);
    for u=1:Ntry        
        imPts = cat(3, pLeft, pRight);
        imPts(1,:,1)=imPts(1,:,1)+randn(1,length(imPts(1,:,1)))*sigma1(v);
        imPts(1,:,2)=imPts(1,:,2)+randn(1,length(imPts(1,:,2)))*sigma1(v);
        imPts(2,:,2)=imPts(2,:,2)+randn(1,length(imPts(2,:,2)))*sigma1(v);
        imPts(1,:,2)=imPts(1,:,2)+randn(1,length(imPts(1,:,2)))*sigma1(v);
        nPts = size(imPts,2);
        if size(imPts,1) == 2
          imPts = cat(1, imPts, ones(1,nPts,2));
        end 
        idTest = randperm(nPts);
        nTest = min(8, nPts);
        idTest = idTest(1:nTest);
        [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);
        [F_1 Sa_1 Sf_1] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),0);

%     %% RANSAC for F
%         seeds = {};
%         sigma = 2.0; rho = 2;
%         for kTrial = 1: nTrial
%           %% Test out F matrix on a random sample of 8 points
%           idTest = randperm(nPts);
%           nTest = min(8, nPts);
%           idTest = idTest(1:nTest);
% 
%           %% Solve for F matrix on the random sample
%           [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);
% 
%           %% Compute perpendicular error of all points to epipolar lines
%           perpErrL = zeros(1,nPts);
%           for k = 1:nPts
%             lk = imPts(:,k,2)' * F';
%             perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
%           end
% 
%           %% Detect inliers
%           idInlier = abs(perpErrL) < rho*sigma;
% 
%           %% Count inliers
%           nInlier = sum(idInlier);
%           if nInlier > 20
%             %% Store sets of sampled points with at least 20 inliers
%             seed = struct;
%             seed.id = idTest;
%             seed.idInlier = idInlier;
%             seed.nInlier = nInlier;
%             seed.F = F;
% 
%             kSeed = length(seeds)+1;
%             seeds{kSeed} = seed;
%           end
%         end 
%         %% Done RANSAC trials
% 
%         %% Extract best solution
%         nInliers = zeros(1, length(seeds));
%         for ks = 1:length(seeds)
%           nInliers(ks) = seeds{ks}.nInlier;
%         end 
%         [nM ks] = max(nInliers);
%         nInliers(ks);
% 
%         %%  Refine estimate of F using all inliers.
%         F = seeds{ks}.F;
%         idInlier = seeds{ks}.idInlier;
% 
%         idInlierOld = idInlier;
%         sum(idInlier);
%         %% Do at most 10 iterations attempting to entrain as many points as possible.
%         for kIt = 1: 10
%           %% Fit F using all current inliers
%           [F Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
%           %% Compute perpendicular error to epipolar lines
%           perpErrL = zeros(1,nPts);
%           for k = 1:nPts
%             lk = imPts(:,k,2)' * F';
%             perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
%           end
%           idInlier = abs(perpErrL) < rho*sigma;
%           nInlier = sum(idInlier);
%           %% If we have the same set of inliers as the previous iteration then stop.
%           if all(idInlier == idInlierOld)
%             break;
%           end
%           idInlierOld = idInlier;
%         end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% RANSAC for F
%         seeds = {};
%         sigma = 2.0; rho = 2;
%         for kTrial = 1: nTrial
%           %% Test out F matrix on a random sample of 8 points
%           idTest = randperm(nPts);
%           nTest = min(8, nPts);
%           idTest = idTest(1:nTest);
% 
%           %% Solve for F matrix on the random sample
%           [F_1 Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),0);
% 
%           %% Compute perpendicular error of all points to epipolar lines
%           perpErrL = zeros(1,nPts);
%           for k = 1:nPts
%             lk = imPts(:,k,2)' * F_1';
%             perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
%           end
% 
%           %% Detect inliers
%           idInlier = abs(perpErrL) < rho*sigma;
% 
%           %% Count inliers
%           nInlier = sum(idInlier);
%           if nInlier > 20
%             %% Store sets of sampled points with at least 20 inliers
%             seed = struct;
%             seed.id = idTest;
%             seed.idInlier = idInlier;
%             seed.nInlier = nInlier;
%             seed.F_1 = F_1;
% 
%             kSeed = length(seeds)+1;
%             seeds{kSeed} = seed;
%           end
%         end 
%         %% Done RANSAC trials
% 
%         %% Extract best solution
%         nInliers = zeros(1, length(seeds));
%         for ks = 1:length(seeds)
%           nInliers(ks) = seeds{ks}.nInlier;
%         end 
%         [nM ks] = max(nInliers);
%         nInliers(ks);
% 
%         %%  Refine estimate of F using all inliers.
%         F_1 = seeds{ks}.F_1;
%         idInlier = seeds{ks}.idInlier;
% 
%         idInlierOld = idInlier;
%         sum(idInlier);
%         %% Do at most 10 iterations attempting to entrain as many points as possible.
%         for kIt = 1: 10
%           %% Fit F using all current inliers
%           [F_1 Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),0);
%           %% Compute perpendicular error to epipolar lines
%           perpErrL = zeros(1,nPts);
%           for k = 1:nPts
%             lk = imPts(:,k,2)' * F_1';
%             perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
%           end
%           idInlier = abs(perpErrL) < rho*sigma;
%           nInlier = sum(idInlier);
%           %% If we have the same set of inliers as the previous iteration then stop.
%           if all(idInlier == idInlierOld)
%             break;
%           end
%           idInlierOld = idInlier;
%         end
        x_vec=linspace(-150,150,30);
        y_vec=linspace(-100,100,30);
        nCol=length(x_vec)^2;
        col = hsv(length(x_vec)^2);
        [xx,yy]=ndgrid(x_vec,y_vec);
        gridpnt=[xx(:) yy(:) ones(length(xx(:)),1)];
        cropBox =[-150  -100   150   100];
        perpErr=[];
        for k= 1:length(gridpnt)
          lk_0 = gridpnt(k,:) * F_0';
          lk = gridpnt(k,:) * F';
          epk = cropLineInBox(lk(1:2), lk(3), cropBox);
          if ~sum(sum(isnan(epk)))
              err=max((lk_0 * [epk(1,:) 1]'/norm(lk_0 (1:2))),(lk_0  * [epk(2,:) 1]'/norm(lk_0 (1:2))));
              perpErr = [perpErr err];
          end
        end
        maxperpErr(u)=max(abs(perpErr));
        perpErr=[];
        for k= 1:length(gridpnt)
          lk_0 = gridpnt(k,:) * F_0';
          lk = gridpnt(k,:) * F_1';
          epk = cropLineInBox(lk(1:2), lk(3), cropBox);
          if ~sum(sum(isnan(epk)))
              err=max((lk_0 * [epk(1,:) 1]'/norm(lk_0 (1:2))),(lk_0  * [epk(2,:) 1]'/norm(lk_0 (1:2))));
              perpErr = [perpErr err];
          end
        end
        maxperpErr1(u)=max(abs(perpErr));
    end
    err_sig(v)=median(abs(maxperpErr));
    err_sig1(v)=median(abs(maxperpErr1));
end

figure(1)
loglog(sigma1,err_sig,'r','LineWidth',2)
hold on
loglog(sigma1,err_sig1,'b','LineWidth',2)
xlabel('\sigma')
ylabel('$d_{\max}$','Interpreter','latex')
axis([10^-5 10^1.6 10^-4 10^4])
ax=gca;
ax.XTick=[ 10^-4 10^-3 10^-2 10^-1 10^0 10^1];
legend('w/ Hartley Norm','w/o Hartley Norm','location','southeast')

figure(2)
plot(sigma1,err_sig,'r','LineWidth',2)
hold on
plot(sigma1,err_sig1,'b','LineWidth',2)
xlabel('\sigma')
ylabel('$d_{\max}$','Interpreter','latex')
legend('w/ Hartley Norm','w/o Hartley Norm','location','southeast')

