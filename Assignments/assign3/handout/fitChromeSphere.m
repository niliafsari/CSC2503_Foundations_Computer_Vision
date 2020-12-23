function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
  mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
  mask = mask(:,:,1) / 255.0;
  mask=round(mask);
  mask=mask';
  ind=find(mask==1);
  [X,Y]=ind2sub(size(mask),ind);
  r=0.25*((max(X)-min(X))+(max(Y)-min(Y)));
  x_c=0.5*(max(X)+min(X));
  y_c=0.5*(max(Y)+min(Y));
  %   min_x=find(mask==1)
  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,n) = im(:,:,1)';           % red channel
  end
  siz=size(imData);
  idx=find(imData==255);
  [xx,yy,z]=ind2sub(siz,idx);
  for n=1:nDir
      x(n)=mean(xx(z==n));
      y(n)=mean(yy(z==n));
  end
  n=[((x-x_c)/r)',((y-y_c)/r)',-(sqrt(r^2-((x-x_c).^2+(x-x_c).^2))/r)'];
  n=n./repmat(sqrt(sum(n.^2,2)),[1 3]); 
  m=repmat([0,0,-1],[nDir,1]);
  L=2*(repmat(dot(n,m,2),[1,3])).*n-m;
  if chatty
    figure
    theta = 0:0.1:2*pi;
    theta = [theta 2*pi];
    plot(cos(theta), sin(theta), 'k');
    hold on
    quiver(((x-x_c)/r)',((y-y_c)/r)',n(:,1),n(:,2),0,'r','lineWidth',2.5);
    hold on
    quiver(((x-x_c)/r)',((y-y_c)/r)',L(:,1),L(:,2),0,'b'); 
    legend('normal','light dir')
    axis equal 
    axis ij
    pause(1);
  end 
  L=L';
  return;

