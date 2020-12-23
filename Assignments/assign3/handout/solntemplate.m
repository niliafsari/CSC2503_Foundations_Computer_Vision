%%
clear 
close all

% Dependencies: 
%   Matlab iseToolbox and utvisToolbox
%   Neet to set the path to include these toolboxes.
%   This is done by running the startup.m file (see below).
if isunix
  % You need to change the following two directories
  % if you are running under unix or linux.
  if ~exist('matlabVisRoot', 'var')
    %run('/home/nilou/Dropbox/Courses/matlab/startup.m');
    run '~jepson/pub/matlab/startup.m';
  end
  #codeDir = '/home/nilou/Dropbox/Courses/Vision/Assignemnts/handout3/';
  codeDir = [pwd, '/'];
  dataDir = [codeDir, 'ppmImages/'];
  chromeDir = [dataDir, 'chrome/'];
else
  % You need to change the following two directories and
  % make sure you have run the startup.m file that
  % came with the utvisToolbox.
  if ~exist('matlabVisRoot', 'var')
   run('C:/Matlab/startup.m');
  end
  codeDir = [pwd, '/'];
  dataDir = [codeDir, 'ppmImages\'];
  chromeDir = [dataDir, 'chrome\'];
end
cd(codeDir)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial coordinates:
% We'll assume a right handed coordinate frame with
% X to the right, Y down, Z in the direction we are looking.
% We assume orthographic imaging, with the camera coords 
% aligned with world coordinates.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
dirMethod  = 1;  % 0 -- Use default light source directions.
                 % 1 -- Use chrome images to estimate lightsource
                 % directions
                 
% The number of different chrome images we have:
nDirChrome = 8;

% The number of different shading images we have (>= nDirChrome)
% the first nDirChrome of which have the same light source
% directions as in the corresponding chrome images.
nDir = 12;

chattyChrome = false;  % show intermediate results in chrome images.
chatty = false;  % Show intermediate results of normal and surface fitting.

% Clear figure to be used for light source directions,
% for which we will superimpose results from all image sets.
figure(2); clf;

% Loop over the test image sets to be used
for useImageSet = 3 % or just choose one image set, e.g. 3
  close all
  switch useImageSet
   case 1
    name = 'gray';
   case 2
    name = 'buddha';
   case 3
    name = 'cat';
   case 4
    name = 'owl';
   case 5
    name = 'horse';
   case 6
    name = 'rock'; 
   otherwise
    fprintf(2, 'Invalid choice of image set, # %d\n', useImageSet);
  end

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 1:  Estimate Light Source Directions using chrome 
  %          sphere images.
  [Lchrome] = getLightDir(dirMethod, chromeDir, nDirChrome, chattyChrome);
  
  % Sanity check
  nrm = sqrt(sum(Lchrome.^2,1));
  if any(abs(nrm - 1) > 1.0e-6)
    fprintf(2, 'Error: Lchrome are not unit vectors\n');
  end
  
  % Plot recovered directions
  theta = 0:0.1:2*pi;
  theta = [theta 2*pi];
  figure(2);
  plot(cos(theta), sin(theta), 'k');
  hold on;
  hLS(1) = plot(Lchrome(1,:), Lchrome(2,:), '*r');
  for k = 1:nDirChrome
    text(Lchrome(1,k), Lchrome(2,k), sprintf(' %d', k-1));
  end
  axis ij; axis equal;
  title('Orthographic image of light source directions.','fontSize',12);
  xlabel('x'); ylabel('y');
  set(gcf, 'Position', [100 100 450 450],'color','w');
%   name=['plots/Ldir.pdf'];
%   export_fig(name, '-pdf');
  if chatty
    fprintf(2, 'Press any key to continue ... ');
    pause; 
  end

  % END OF PART 1.
  % Uncomment/comment the following continue statement.  This
  % allows you to run your partially completed code on
  % each of the examples.
  %continue;
  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read in images of the object

  imageDir = [dataDir, name, '/'];

  % Read mask image and binarize
  mask = ppmRead([imageDir,name,'.mask.ppm']);
  mask = mask(:,:,1) / 255.0;
  mask = mask > 0.5;

  % Get image size and number of pixels.
  imsize = [size(mask,1), size(mask,2)];
  numPixels = prod(imsize);

  % Vectorize mask
  mask = mask(:);

  % We will switch between storing images as 2D arrays of size
  % imsize(1) x imsize(2) and vectorizing them as long numPixels x 1
  % vectors.  Here, by default, we store them in the vector form,
  % since this makes most of the operations we need to do easier
  % (other than display).  When we integrate the normals to get z
  % it will be convenient to reshape the normal image to be a
  % imsize(1) x imsize(2) x 3 array.

  % Read in images, store gray-scale vectorized images.
  imData = zeros(numPixels, nDir);
  for n=1:nDir
    fname = [imageDir,name,'.',num2str(n-1),'.ppm'];
    RGBim = ppmRead(fname);
    if chatty
      figure(1); clf;
      image(RGBim/255); 
      axis equal;
      axis off;
      title(sprintf('Data image %d',n));
      pause(0.5);
    end
    % Calculate a grayscale image
    imGray = sum(RGBim,3)/3;
    imData(:, n) = imGray(:);
  end


  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 2: Fit surface normals and the albedo for pixels within 
  %         the object mask using the images for which corresponding 
  %         images of the chrome sphere (and hence light source 
  %         directionss) are available.

  imDataCrop = imData(mask,1:nDirChrome);

  [nCrop, albedoCrop] = fitReflectance(imDataCrop, Lchrome);


  % Unpack the normals and albedos estimated from within the mask
  % into imsize sized images, and display.

  n = zeros([numPixels, 3]);
  n(mask,:) = nCrop;
  albedoGray = zeros(numPixels,1);
  albedoGray(mask) = albedoCrop;

  % Display gray albedo
  figure(3); clf;
  showIm(reshape(albedoGray,imsize) );
  title('Recovered albedo (gray)');
  pause(1);
%     set(gcf, 'Position', [100 100 700 550],'color','w');
%     name1=['plots/' name 'alb.png'];
%     export_fig(name1, '-png');
  % Display each component of the normal as a separate image.
  n = reshape(n, [imsize, 3]);
  figure(4)
  subplot(2,2,1);
  showIm(n(:,:,1));
  title('Surface Normal (nx)');
  subplot(2,2,2);
  showIm(n(:,:,2));
  title('Surface Normal (ny)');
  subplot(2,2,3);
  showIm(n(:,:,3));
  title('Surface Normal (nz)');
  n = reshape(n, [numPixels, 3]);
  pause(1);
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=['plots/'  name 'norm.png'];
%   export_fig(name1, '-png');
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute the reconstruction error and display
   L = Lchrome;
  rmsErr = zeros(numPixels, 1);
  for k = 1:nDirChrome
    nDotL = n * L(:,k);
    rec = nDotL .* albedoGray;
    err = (rec -  imData(:,k)).*mask;
    f=figure(5); clf;
    subplot(2,2,1); showIm(reshape(imData(:,k), imsize));
    title(sprintf('Image %d', k));
    subplot(2,2,2); showIm(reshape(rec, imsize));
    title('Reconstruction');
    subplot(2,2,3); showIm(reshape(err, imsize));
    title('Error');
    subplot(2,2,4); showIm(reshape(double(nDotL <0), imsize));
    title('n dot L < 0');
    rmsErr = rmsErr + err.^2;
%    set(f, 'Position', [100 100 750 500],'color','w');
    pause(1);
%     name1=[ 'plots/'  name  'rec_' mat2str(k) '.png'];
%     export_fig(name1, '-png');
    
  end
  rmsErr = sqrt(rmsErr/nDirChrome);
  figure(6); clf; showIm(reshape(rmsErr, imsize));
  title(sprintf('RMS error (total: %f)', sqrt(sum(rmsErr.^2))));
  pause(1);
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots/'  name 'err.png'];
%   export_fig(name1, '-png');
  if chatty
    fprintf(2, 'Press any key to continue ... ');
    pause; 
  end

  % END OF PART 2.
  % Uncomment/comment the following continue statement.  This
  % allows you to run your partially completed code on
  % each of the examples.
  %continue;
  

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 3: Given normals and light source directions, 
  %         calculate the albedo in each color channel
  RGBim = zeros(numPixels,3, nDirChrome);
  RGBimcrop = zeros(sum(mask),3, nDirChrome);
  albedo = zeros(numPixels, 3);
  for j=1:nDirChrome
    fname = [imageDir,name,'.',num2str(j-1),'.ppm'];
    temp = ppmRead(fname);
    RGBim(:,:,j)=reshape(temp, numPixels,3);
  end 
  for i=1:3
    RGBimcrop(:,i,1:nDirChrome) = RGBim(mask,i,1:1:nDirChrome);
    im=reshape(RGBimcrop(:,i,:),sum(mask),nDirChrome);
    g=nCrop*L;
    albedo(mask,i)= dot(g,im,2)./dot(g,g,2);
  end
  

  % YOU NEED TO ADD CODE HERE FOR PART 4

 
  % Clip albedo to range [0, 255] 
  albedo = max(albedo,0);
  albedo = min(albedo,255);
  figure(3); clf; 
  image(reshape(albedo/255, [imsize, 3]));
  title('Recovered Albedo (RGB)');
  axis equal; axis off;
  pause(1);
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots3/'  name 'rec_alb.png'];
%   export_fig(name1, '-png');
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Show images for synthetic light sources
  t = 0:0.15:2*pi;
  r = 0.5;
  Lsyn = zeros(3, length(t));
  Lsyn(1,:) = r * cos(t);
  Lsyn(2,:) = r *sin(t);
  Lsyn(3,:) = -sqrt(1 - sum(Lsyn(1:2,:).^2,1));

  figure(2); hold on;
  hLS(2) = plot(Lsyn(1,:), Lsyn(2,:), '-*g');

  for k = 1:size(Lsyn,2)
    im = albedo .* repmat(n * Lsyn(:,k), 1,3);
    im = max(im,0);
    im = min(im, 255);
    figure(1); clf; 
    image(reshape(im/255, [imsize, 3]));
    axis equal; axis off;
    title('Synthetically Shaded Image');
%     if ~mod(k,10)
%         set(gcf, 'Position', [100 100 700 550],'color','w');
%         name1=[ 'plots3/'  name '_' mat2str(k) '_syn.png'];
%         export_fig(name1, '-png');
%         figure(2); hold on;
%         k
%         hLS(3) = plot(Lsyn(1,k), Lsyn(2,k), '-*m');
%     end
    pause(0.5);
  end
%   
  figure(2); hold on;
%  legend(hLS(1:3), {'Chrome', 'Synthetic','shown in report'});
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots3/syn.png'];
%   export_fig(name1, '-png');

  if chatty
    fprintf(2, 'Press any key to continue ... ');
    pause; 
  end

  % END OF PART 3.
  % Uncomment/comment the following continue statement.  This
  % allows you to run your partially completed code on
  % each of the examples.
  %continue;


  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 4: Estimate depth from normals
  maskDepth = mask & (sum(n.^2,2) > 0);
  [depth] = getDepthFromNormals(reshape(n, [imsize, 3]), ...
                                reshape(maskDepth,imsize));

  % Display depth map as image and as surface mesh.
  % To put the mesh in the correct perspective, invert depth and translate
  % (because smaller depth is closer to the camera, thus 'higher'... )
  figure(5); 
  displayImage(depth);
  title('Estimated Object Depth');
%      set(gcf, 'Position', [100 100 700 550],'color','w');
%     name1=[ 'plots4/'  name '_depth.png'];
%     export_fig(name1, '-png');
  if 1
    % On some Unix/Linux machines, the following
    % code can crash.  I believe the problem is with meshz().
    s = depth + 1.01*max(depth(:))*reshape(~maskDepth,imsize);
    s = max(s(:))-s;
    s = fliplr(s);
    figure(6);
    clf;
    meshz(s);
    colormap(jet(256));
    title('Estimated Depth as Mesh');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
%      set(gcf, 'Position', [100 100 700 550],'color','w');
%     name1=[ 'plots4/'  name '_depthmesh.png'];
%     export_fig(name1, '-png');
    pause(1);
  end

  if chatty
    fprintf(2, 'Press any key to continue ... ');
    pause; 
  end

  % END OF PART 4.
  % Uncomment the following continue statement.  This
  % allows you to run your partially completed code on
  % each of the examples.
  %continue;


  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 5:  Estimate light source directions for four 
  % BONUS    additional images
  g = nCrop .* repmat(albedoCrop, 1, 3);
  Lrec = zeros(3,nDir-nDirChrome);

  Lrec = inv(g'*g)*g'*imData(mask,(nDirChrome+1):nDir);
  

  LrecMag = sqrt(sum(Lrec.^2,1));
  figure(5); clf;
  plot(LrecMag, '-*b'); hold on;
  plot(ones(length(LrecMag), 1), '-k');
  title('Recovered Additional Light Source Magnitudes');
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots5/'  name '_recadd.png'];
%   export_fig(name1, '-png');
  LrecDir = Lrec ./ repmat(LrecMag, 3, 1);
  figure(2); hold on;
  hLS(2) = plot(LrecDir(1,:), LrecDir(2,:), '*b');
  for k = 1:size(Lrec,2)
    text(LrecDir(1,k), LrecDir(2,k), ...
         sprintf(' %d:%d', k+nDirChrome, nDirChrome+useImageSet));
  end
  legend(hLS(1:2), {'Chrome', 'Recovered'});
%  set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots5/'  name '_recL.png'];
%   export_fig(name1, '-png');
  if chatty
    fprintf(2, 'Press any key to continue ... ');
    pause;
  end

  % Concatenate chrome and recovered light source directions.
  L = [Lchrome Lrec];

  % Compute the reconstruction error and display
  rmsErr = zeros(numPixels, 1);
  for k = 1:nDir
    nDotL = n * L(:,k);
    rec = nDotL .* albedoGray;
    err = (rec -  imData(:,k)).*mask;
    figure(5); clf;
    subplot(2,2,1); showIm(reshape(imData(:,k), imsize));
    title(sprintf('Image %d', k));
    subplot(2,2,2); showIm(reshape(rec, imsize));
    title('Reconstruction');
    subplot(2,2,3); showIm(reshape(err, imsize));
    title('Error');
    subplot(2,2,4); showIm(reshape(double(nDotL <0), imsize));
    title('n dot L < 0');
    rmsErr = rmsErr + err.^2;
    pause(1);
  end
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots5/'  name '_part5err.png'];
%   export_fig(name1, '-png');
  rmsErr = sqrt(rmsErr/nDir);
  figure(6); clf; showIm(reshape(rmsErr, imsize));
  title(sprintf('RMS error (total: %f)', sqrt(sum(rmsErr.^2))));
  pause(2);
%   set(gcf, 'Position', [100 100 700 550],'color','w');
%   name1=[ 'plots5/'  name '_part5.png'];
%   export_fig(name1, '-png');
  %%%%%%%%%%%%%%%%%%%%%%%%%%
end % loop over different objects at beginning of script
