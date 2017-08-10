% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function val = supermodularZeroSingleton(func,x);
val = func(x);
for i=1:length(x)
  val = val - func(x(i));
end

end
