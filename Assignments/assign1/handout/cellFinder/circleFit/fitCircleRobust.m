function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.


% initialization
siz=size(pts);
i=1;
eps=20000;
N=siz(1);
P=zeros(3,1);
P(1:2,1)=-2*initx0;
P(3,1)=initx0(1)^2+initx0(2)^2-initr^2;
U=pts.^2;
V=[pts ones(N,1)];
D=[U V];
niter=500;


while (eps>0.01 && i<niter)
    e=D*[ones(2,1); P];
    % weights
    W=((2*sigmaGM^2*e)./(sigmaGM^2+e.^2))./e;
    x0=-P(1:2,1)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %finding leverage points
    threshold = prctile(W,85);
    high_pts=pts(W>threshold,:);
    W((W/max(W))< 0.016,:)=0;
    if ~isempty(high_pts)
        [theta,~]=cart2pol(x0(1)-high_pts(:,1),x0(2)-high_pts(:,2));
        theta(theta<0)=2*pi+theta(theta<0);
        [theta,id]=sort(theta,'ascend');
        high_pts=high_pts(id,:);
        all_ids=1:length(high_pts);
        theta_diff=zeros(size(theta));
        theta_diff(2:length(theta),1)=theta(2:length(theta))-theta(1:length(theta)-1);
        theta_diff(1,1)=theta(1,1)-theta(length(theta),1);
        theta_diff(theta_diff<0)=2*pi+theta_diff(theta_diff<0);
        theta_diff(theta_diff>2*pi)=theta_diff(theta_diff>2*pi)-2*pi;
        i_end=find((theta_diff-median(theta_diff))>(2.5*std(theta_diff)));
        %is there is any gap:
        if ~isempty(i_end)
            i_start=i_end-1;
            i_start(i_start==0)=length(theta_diff);
            theta=theta-theta(i_start(1));
            theta(theta<=0)=theta(theta<=0)+2*pi;        
            gaps=[i_start i_end];
            [~,b]=sort(gaps(:,1),'ascend');
            gaps=gaps(b,:);
            siz_g=size(gaps);
            num_gaps=siz_g(1);
            if num_gaps>1
               for u=1:num_gaps
                   if u==num_gaps
                       ids{u}= find(theta>=theta(gaps(u,2)) & theta<=theta(gaps(1,1)));
                   else
                       ids{u}= find(theta>=theta(gaps(u,2)) & theta<=theta(gaps(u+1,1)));
                   end
                   if u==1
                       siz_u=length(ids{u});
                       rm_u=u;
                   end
                   if siz_u<length(ids{u})
                       siz_u=length(ids{u});
                       rm_u=u;
                   end                  
               end
               toremove=setdiff(all_ids,ids{rm_u});
               % setting the weights of leverage points to zero
               for k=1:length(toremove)
                   [~,ii]=min((high_pts(toremove(k),1)-pts(:,1)).^2+(high_pts(toremove(k),2)-pts(:,2)).^2);
                   W(ii)=0;
               end
            end        
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diagonal matrix of weights
    W_diag=diag(W);
    % solve the linear equation and obtain new parameters
    P_new=linsolve(V'*W_diag*V,-V'*W_diag'*U*ones(2,1));
    %updating parameters and epsilon
    eps=norm(P_new-P);
    P=P_new;
    x0=-P(1:2,1)/2;
    i=i+1;
end
%ouputs
r=sqrt(norm(x0)^2-P(3,1));
w= W;
maxW=max(w);

