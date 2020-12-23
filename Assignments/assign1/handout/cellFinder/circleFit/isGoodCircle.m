function [goodCircle] = isGoodCircle(x0, r, w,...
                                     circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.
  
  if ~ nFound
     goodCircle=true; %for the first circle always true
  else
      for i=1:nFound
          d=sqrt((x0(1)-circleEstimates(1,i))^2+(x0(2)-circleEstimates(2,i))^2);
          if (r+circleEstimates(3,i))<=d %circles are not overlapping
              continue
          else
             if d < abs(circleEstimates(3,i) -r) % one circle is contained in another
                goodCircle=false;
                return
             else % overlapping cases
                A=area_intersect_circle_analytical([x0(1) x0(2) r; circleEstimates(1,i) circleEstimates(2,i) circleEstimates(3,i)]);
                % we don't allow more that 38% overlap for circles larger
                % than 20% of a previously found circle
                if (((A(2,1)/(pi*r^2))>0.38) || ((A(2,1)/(pi*circleEstimates(3,i)^2))>0.38)) && ((A(1,1)/A(2,2))>0.2)
                    goodCircle=false;
                    return 
                end
             end
          end
      end
      % compute the sum of intersection with all previous circles
      A=area_intersect_circle_analytical([x0(1) x0(2) r; circleEstimates(:,1:nFound)']);
      A=sum(A(1,:))-A(1,1);
      % reject the proposed circle if sum if larger than 40% of its surface
      if ((A/(pi*r^2))>0.4)
        goodCircle=false;
        return 
      end
      %in all other cases return true!
      goodCircle=true;      
  end
  
  % YOU FINISH THIS