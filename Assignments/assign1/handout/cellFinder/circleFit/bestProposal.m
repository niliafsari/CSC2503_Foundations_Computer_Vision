function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals
  
  % YOU COMPLETE THIS
  
  %initialization
  siz_data=size(p);
  siz_cir=size(circles);
  e=zeros(siz_cir(1),siz_data(1));
  gm=zeros(siz_cir(1),siz_data(1));
  similarity=zeros(siz_cir(1),siz_data(1));
  %calculating gm and similarity measures
  for i=1:siz_cir(1)
      for j=1:siz_data(1)
          e(i,j)=[p(j,1)^2 p(j,2)^2 p(j,1) p(j,2) 1]*[1 ; 1; -2*circles(i,1) ...
              ;-2*circles(i,2); circles(i,1)^2+circles(i,2)^2-circles(i,3)^2];
          gm(i,j)=e(i,j)^2/(e(i,j)^2+sigmaGM^2);
          similarity(i,j)=dot([circles(i,1);circles(i,2)]-[p(j,1); p(j,2)],normals(j,:)')/norm([p(j,1);p(j,2)]-[circles(i,1);circles(i,2)]);
      end
  end
  % normalization
  normalized_gm=(sigmaGM^2*gm-1)/(sigmaGM^2-1);
  normalized_similarity=(-similarity+1)/2;
  normalized_gm=mean(normalized_gm,2);
  normalized_similarity=mean(normalized_similarity,2);
  normalized_gm=(normalized_gm-mean(normalized_gm))/std(normalized_gm);
  normalized_similarity=(normalized_similarity-mean(normalized_similarity))/std(normalized_similarity);
  % computing rank measure
  idx=1:siz_cir(1);
  proposed_rank=(idx/siz_cir(1));
  proposed_rank=(proposed_rank-mean(proposed_rank))/std(proposed_rank);
  % combining all measure to come up with one cost function
  cost=0.5*mean(normalized_gm,2)+1.5*mean(normalized_similarity,2)+0.5*proposed_rank';
  % pick the minimum cost
  [~,index]=sort(cost,'ascend');  
  if siz_cir(1)>=4
    circle=circles(index(randi([1,4])),:); 
  else
    circle=circles(index(randi([1,siz_cir(1)])),:); 
  end