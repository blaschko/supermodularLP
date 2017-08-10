% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918

function [ind] = randomSelection (func, d, budget)

ind = [ones(budget,1);zeros(d-budget,1)];
ind = find(ind(randperm(length(ind))));

end
