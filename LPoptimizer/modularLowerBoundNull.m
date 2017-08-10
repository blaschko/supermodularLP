function m = modularLowerBoundNull(g,p);
% Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% copyright 2016
%
% Implements m g from Definition 8 of
% Blaschko, M. B.: Slack and Margin Rescaling as Convex Extensions of
% Supermodular Functions. Version 1 arXiv:1606.05918, 2016.
%
%


for i=1:p
    A = i;
    m(i) = g(A);
end

end
