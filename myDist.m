% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918

function [D] = myDist (X, Y)
  if(nargin<2)
    Y=X;
  end 
  D = sqrt(sum(X.*X,2)*ones(1,size(Y,1)) + ones(size(X,1),1)*sum(Y.*Y,2)' - 2*X*Y');
end
