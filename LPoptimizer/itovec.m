% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function curvec = itovec(i,p)
  assert(i<2^p) % check for semantically valid input
  assert(i>=0)
  assert(sum(itovec2((2^p)-1,p))==p) %check for overflow
  curvec = itovec2(i,p);
end
    
function curvec = itovec2(i,p)
curvec = zeros(p,1);
for j=1:p
    curvec(j) = mod(i,2);
    i = floor(i/2);
end
end
