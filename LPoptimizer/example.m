% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function [x,value,LPvalue,func,state,bruteForceMinimum] = example(slackmargin,b,v,budget);

if(exist('budget','var'))
    sparm.budget = budget;
end

if(nargin<2)
    b = 2;
end
sparm.p = 2*b;
if(nargin<3)
    v = zeros(sparm.p,1);
end
sparm.func = @(A)(theFunc(A,b,v));
gamma = findMaxGammaMarginRescaling(sparm.func,sparm.p);
sparm.func = @(A)(gamma*theFunc(A,b,v));
if(nargin==0)
    sparm.formulationType = 'slack';
else
    sparm.formulationType = slackmargin;
end

assert(isSupermodular(sparm.func,sparm.p));

[x,sparm,state] = LPcuttingPlane(sparm);
value = sparm.func(find(round(x)));
LPvalue = sparm.bestobjective;

func = sparm.func;

if(nargout>=6)
    A = setfnMaxBruteForce(@(A)(-sparm.func(A)),sparm.p)
    bruteForceMinimum = A;
end

end


function val = theFunc(A,b,v)
  if(nargin<2)
      b = 3;
  end
  if(nargin<3)
      v = zeros(b*2,1);
  end

  sizeA = length(unique(A));

  val = (((sizeA-b)^2)-b^2)/(b*b-(b-1)*(b-1));
  xA = zeros(size(v));
  xA(A) = 1;
  val = val + dot(v,xA);
end



