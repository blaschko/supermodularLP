% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918

function [func,D,coverageInd] = coverageFunction(data,radiusThresh)

D2 = myDist(double(data));

if(nargin<2)
  radiusThresh = quantile(D2(:),0.1);
end

D2 = (D2<radiusThresh);

if(nargout>=3)
  % remove duplicates
  isDuplicate = myDist(double(D2))==0;

  for i=1:size(isDuplicate,2)
    sameInd = find(isDuplicate(:,i));
    if(min(sameInd)<i)
      coverageInd(i) = 0;
    else
       coverageInd(i) = 1;
    end
  end

  coverageInd = find(coverageInd);

  D = D2(coverageInd,:);
else
  D = D2;
end

func = @(A)(-sum(sum(D(A,:),1)>0,2));

end
