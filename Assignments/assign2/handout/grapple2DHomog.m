%% File: grappleFmatrix
%% A2 2017 handout code
%% Uses RANSAC to estimate F matrix from corresponding points.
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
   cd ~jepson/pub/matlab   %% CHANGE THIS
   startup;
   cd(dir);
end

reconRoot = '.';  
addpath([reconRoot '/data/wadham']);
addpath([reconRoot '/utils']);


% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  % Number of ransac trials to try.

%% Wadham left image: use  wadham/001.jpg
imPath = 'data/wadham/'; fnameLeft = '001'; 
im = imread([imPath fnameLeft],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imLeft = imDwn;

%% Read correspondence data
load data/wadham/corrPnts2
%% Wadham right image data/wadham002-5.jpg use for corrPnts2-5 respectively
fnameRight = '002';
im = imread([imPath fnameRight],'jpg');
im = double(im(:,:,2));
imDwn = blurDn(im,1)/2;
imRight = imDwn;

clear imPts;
imPts = cat(3,  im_pos2',im_pos1');
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end

%% RANSAC for F
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  %% Test out F matrix on a random sample of 8 points
  idTest = randperm(nPts);
  nTest = min(4, nPts);
  idTest = idTest(1:nTest);

  %% Solve for F matrix on the random sample
  [H, ~] = linEstH(imPts(:,idTest,1), imPts(:,idTest,2),1);
  
  %% Compute error of all points to the correspondences
  errL = zeros(1,nPts);
  for k = 1:nPts
    lpnt = H*imPts(:,k,2);
    lpnt=lpnt/lpnt(3);
    errL(k) = norm(lpnt-imPts(:,k,1));
  end
  
  %% Detect inliers
  idInlier = abs(errL) < rho*sigma;
  
  %% Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    %% Store sets of sampled points with at least 20 inliers
    seed = struct;
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.H = H;
    
    kSeed = length(seeds)+1;
    seeds{kSeed} = seed;
  end
end 
%% Done RANSAC trials

%% Extract best solution
nInliers = zeros(1, length(seeds));
for ks = 1:length(seeds)
  nInliers(ks) = seeds{ks}.nInlier;
end 
[nM ks] = max(nInliers);
nInliers(ks)

%%  Refine estimate of F using all inliers.
H = seeds{ks}.H;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  %% Fit F using all current inliers
  [H,~] = linEstH(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
  
  %% Compute perpendicular error to epipolar lines
  errL = zeros(1,nPts);
  for k = 1:nPts
    lpnt = H*imPts(:,k,2);
    lpnt=lpnt/lpnt(3);
    errL(k) = norm(lpnt-imPts(:,k,1));
  end
  idInlier = abs(errL) < rho*sigma;
  nInlier = sum(idInlier)
  
  %% If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end
  
%%%%%%%%%%%%%%%%%%%%% Plot results
nTest = 64;  %% Number of epipolar lines to plot
nCol = 16;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

%% Random sample the lines to plot
idLines = randperm(nPts);  
idLines = idLines(1:nTest);

%% Show left image
SUPERIMPOSE = TRUE;
hFig = figure(1);
clf;
if SUPERIMPOSE
  image(imRight);
  colormap(gray(256));
end
resizeImageFig(hFig, size(imRight), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');



SUPERIMPOSE = TRUE;
hFig = figure(2);
clf;
if SUPERIMPOSE
  image(homogWarp(imLeft,inv(H)));
  colormap(gray(256));
end
resizeImageFig(hFig, size(imRight), 1); hold on;
set(get(hFig,'CurrentAxes'),'Ydir','reverse');

errL = [];
for k = 1:nPts
  lpnt = H*imPts(:,k,2);
  lpnt=lpnt/lpnt(3);
  errL = [errL norm(lpnt-imPts(:,k,2))];
end
errR = [];
for k = 1:nPts
  rpnt = inv(H)*imPts(:,k,1);
  rpnt=rpnt/rpnt(3);
  errR = [errR  norm(rpnt-imPts(:,k,2))];
end

%% Plot a histogram of the perpendicular distances
err = [errL errR];
err = min(err, 10);
err = max(err, -10);
[n b] = histo(err, 64);
figure(3); clf;
plot(b,n);
title('Distance Differences for Homography');
xlabel('Error in pixels');
ylabel('Frequency');

%% Count inliers
idL = abs(errL)< rho*sigma;
idR = abs(errR) < rho*sigma;
idInlier = idL & idR;
sum(idInlier)
sum(idInlier)/nPts

