function [depth] = getDepthFromNormals(n, mask)
  % [depth] = getDepthFromNormals(n, mask)
  %
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %
  siz = size(mask);
  M=siz(2);
  N=siz(1);
  size(n)
  n=reshape(n,M*N,3);
  U=sparse(M*N,M*N);
  Q=sparse(M*N,M*N);
  V=zeros(2*M*N,1);
  for i=1:N
      for j=1:M
          if mask(i,j)
              idx=sub2ind(siz,i,j);
              U(idx,idx)=-n(idx,3);
              Q(idx,idx)=-n(idx,3);
              V(idx)=-n(idx,2);
              V(idx+M*N)=-n(idx,1);
              if ((i+1)<=N)
                idx_x=sub2ind(siz,i+1,j);
                U(idx,idx_x)=n(idx,3);   
              end
              if ((j+1)<=M)
                idx_y=sub2ind(siz,i,j+1);
                Q(idx,idx_y)=n(idx,3);   
              end               
          end
      end
  end
  A=[U;Q];
  z = A\V;
  depth = reshape(z, N, M);
  depth(mask==0) = 0;
  % YOU NEED TO COMPLETE THIS.
