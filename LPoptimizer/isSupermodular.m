% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function val = isSupermodular(func,p);

A = submodConstraints(p);

v = zeros(2^p,1);

for i=1:(2^p)
    v(i) = func(find(itovec(i-1,p)));
end

ind = find(A*v<0);
val = (length(ind)==0);

end
