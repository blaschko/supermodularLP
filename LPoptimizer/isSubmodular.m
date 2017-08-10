% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function val = isSubmodular(func,p);

val = isSupermodular(@(x)(-func(x)),p);

end
