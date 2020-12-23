function [circles] = getProposals(normals, p, numGuesses)
  % [circles] = getProposals(normals, p, numGuesses)
  % Attempt to produce up to numGuesses circle proposals from
  % the edgel data p and normals.  For typical data sets
  % we will be able to produce numGuesses proposals.  However,
  % on some datasets, say with only a few edgels, we may not be
  % able to generate any proposals.  In this case size(circles,1)
  % can be zero.
  % Input:
  %  normals - N x 2 edgel normals
  %  p         N x 2 edgel positions
  %  numGuesses - attempt to propose this number of circles.
  % Return:
  %   circles a P x 3 array, each row contains [x0(1) x0(2) r]
  %           with 0 <= P <= numGuesses.  

  % initialize and set the ranges and resolution for hough transform matrix
  distances_origin=sqrt(p(:,1).^2+p(:,2).^2);
  d_max= max(distances_origin);
  d_min= min(distances_origin);
  xc_max=max(p(:,1))+d_max-d_min;
  xc_min=min(p(:,1))-(d_max-d_min);
  yc_max=max(p(:,2))+d_max-d_min;
  yc_min=min(p(:,2))-(d_max-d_min);
  xc_bin=80;
  yc_bin=80;
  r_bin=40;
  r_vec=linspace(0,d_max-d_min,r_bin);
  xc_vec=linspace(xc_min,xc_max,xc_bin);
  yc_vec=linspace(yc_min,yc_max,yc_bin);
  param_space=zeros(r_bin,xc_bin,yc_bin);
  N=length(p);
  sample=N;
  xc_gen=repmat(p(:,1),1,r_bin)+normals(:,1)*r_vec;
  yc_gen=repmat(p(:,2),1,r_bin)+normals(:,2)*r_vec;
  rc_gen=repmat(r_vec,sample,1);
  xc_gen=reshape(xc_gen,[r_bin*sample 1]);
  yc_gen=reshape(yc_gen,[r_bin*sample 1]);
  rc_gen=reshape(rc_gen,[r_bin*sample 1]);
  % parsing the parameter space
  for i=1:r_bin*sample
     [~,i_xc]=min(abs(xc_gen(i)- xc_vec));
     [~,i_yc]=min(abs(yc_gen(i)- yc_vec));
     [~,i_rc]=min(abs(rc_gen(i)- r_vec));
     % cast the vote in the corresponding elements of the parameter matrix
     param_space(i_rc,i_xc,i_yc)=param_space(i_rc,i_xc,i_yc)+1;
  end
  siz=size(param_space);
  param_space=reshape(param_space,[siz(1)*siz(2)*siz(3) 1]);
  % sorting based on the number of votes
  [~,ind]=sort(param_space,'descend');
  I=ind(1:numGuesses);
  [u,j,k]=ind2sub([siz(1) siz(2) siz(3)],I);
  % output based on the most votes
  circles=[ xc_vec(j)' yc_vec(k)' r_vec(u)'];
