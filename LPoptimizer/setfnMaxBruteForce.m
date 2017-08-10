% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function A = setfnMaxBruteForce(g,p)

bestVal = -Inf;
for i=0:2^p-1
  currA = find(itovec(i,p));
  if(g(currA)>bestVal)
    bestVal = g(currA);
    A = currA;
  end
end

end
