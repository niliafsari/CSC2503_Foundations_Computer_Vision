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
  cd ~jepson/pub/matlab   %% CHANGE THIS to your startup directory
   %% CHANGE THIS to your startup directory
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
Ntry=15;
sigma1=logspace(-5,1.65,80);
err_sig=zeros(1,length(sigma1));
for v=1:length(sigma1)
    v
    maxperpErr=zeros(1,Ntry);
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
        maxperpErr(u)=max(perpErr);
    end
    err_sig(v)=median(maxperpErr);
end
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 0.1);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 0.1);
Rright = MextRight(:, 1:3);

sclZ = 0.1;
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
Ntry=15;
err_sig1=zeros(1,length(sigma1));
for v=1:length(sigma1)
    v
    maxperpErr=zeros(1,Ntry);
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
        maxperpErr(u)=max(perpErr);
    end
    err_sig1(v)=median(maxperpErr);
end
figure(1)
loglog(sigma1,err_sig,'r','LineWidth',2)
hold on
loglog(sigma1,err_sig1,'b','LineWidth',2)
xlabel('\sigma')
ylabel('$d_{\max}$','Interpreter','latex')
axis([10^-5 10^1.6 10^-4 10^4])
ax=gca;
ax.XTick=[10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1];
legend('Scz=1','Scz=0.1','location','southeast')

plot(sigma1,err_sig,'r','LineWidth',2)
hold on
plot(sigma1,err_sig1,'b','LineWidth',2)
xlabel('\sigma')
ylabel('$d_{\max}$','Interpreter','latex')
legend('Scz=1','Scz=0.1','location','southeast')
