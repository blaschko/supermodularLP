% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function i = vectoi(curvec)
%warning('not currently checking for overflows - fix')
    j = 2.^int64([0:length(curvec)-1])';
    i = sum(int64(curvec).*j);
end

