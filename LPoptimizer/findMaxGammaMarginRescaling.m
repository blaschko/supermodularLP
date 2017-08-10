function [gamma,A,B] = findMaxGammaMarginRescaling(tildeg,p);
% Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% copyright 2016
% finds maximal value of gamma satisfying Equation (31) of
% Blaschko, M. B.: Slack and Margin Rescaling as Convex Extensions of
% Supermodular Functions. Version 1, arXiv:1606.05918, 2016.
%
% inputs:
%  tildeg -  (interpret as latex \tilde{g}) a supermodular function
%  p - the size of the base set, |V|
%
% returns gamma=1 for constant functions
%

m = modularLowerBoundNull(tildeg,p);

g = @(A)(tildeg(A) - sum(m(A)));

gV = g(1:p);

best = Inf;
i = -1;
for j=1:p
	disp(sprintf(' findMaxGammaMarginRescaling: gVmj, j=%d %s',j, datestr(now)))
	
  gVmj = g([1:j-1 j+1:p]);
  if(gVmj<best)
    i = j;
    best = gVmj;
  end
end

gamma = 1/(gV - best);
  
if(gV==0)
  gamma = 1;
end

B = 1:p;
A = [1:i-1 i+1:p];

end
