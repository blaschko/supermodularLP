% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function [phi,b] = computeJ2constraint(sparm,func,A)
  % computes gradient of bound for linear extension derived from bound 1 in Section 2.2 of R. Iyer and J. Bilmes. Algorithms for approximate minimization of the difference between submodular functions, with applications. In UAI, pages 407–417, 2012.
  %xi >= b - <phi,x>

  fA = func(A);
  Ac = ones(sparm.p,1);
  Ac(A) = 0;
  Ac = find(Ac);

  b = fA;
  phi = zeros(sparm.p,1);
  for i=1:length(Ac)
    phi(Ac(i)) = -func(Ac(i));
  end

  for i=1:length(A)
    phi(A(i)) = -fA + func(A([1:i-1 i+1:end]));
    b = b + phi(A(i));
  end
end
