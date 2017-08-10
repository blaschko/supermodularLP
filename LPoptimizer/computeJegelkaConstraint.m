% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function [phi,b] = computeJegelkaConstraint(sparm,func,A)
  % computes gradient of bound for Equation (52) from Blaschko, M. B.: Slack and Margin Rescaling as Convex Extensions of Supermodular Functions. arXiv:1606.05918 version 1, 2016.
  %xi >= b - <phi,x>
  fV = func(1:sparm.p);
  fA = func(A);
  Ac = ones(sparm.p,1);
  Ac(A) = 0;
  Ac = find(Ac);

  b = func(A);
  phi = zeros(sparm.p,1);
  for i=1:length(Ac)
    phi(Ac(i)) = fA - func(sort([A Ac(i)]));
  end
  for i=1:length(A)
    phi(A(i)) = -fV + func([1:A(i)-1 A(i)+1:sparm.p]);
    b = b + phi(A(i));
  end
end
