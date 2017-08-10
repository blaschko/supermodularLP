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
% Currently implemented very naively
% TODO: determine whether a faster, possibly poly-time algorithm exists

m = modularLowerBoundNull(tildeg,p);

g = @(A)(tildeg(A) - sum(m(A)));

gamma = Inf;

if(((2^p)^2)>1e7)
    A = 1:p;
    gamma = 1/g(A);
    B = zeros(size(A));
    return
end

for i=0:((2^p)-1)
    for j=0:((2^p)-1)
        A = itovec(i,p);
        B = itovec(j,p);
        if(g(find(A))>g(find(B)))
            tmp = A;
            A = B;
            B = tmp;
        end
        if(g(find(B))>g(find(A)))
	  if(gamma>(dot(B,B) - dot(A,B))/(g(find(B))-g(find(A))))
            gamma = min(gamma,(dot(B,B) - dot(A,B))/(g(find(B))-g(find(A))));
            bestB = B;
            bestA = A;
          end
        end
    end
end

if(isinf(gamma))
    gamma = 1;
end

B = bestB;
A = bestA;

end
