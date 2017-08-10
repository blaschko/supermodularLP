% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function ind = coverageNemhauserWolsey(D,budget)
% implements greedy maximization of coverage under a budget.  D is a binary matrix such that D(i,j)=1 iff i covers j.  Nemhauser, G.L., Wolsey, L.A., Fisher, M.L.: An analysis of approximations for maximizing submodular set functions. Mathematical Programming 14, 265–294 (1978)

inds = 1:size(D,1);
ind = [];
covered = zeros(1,size(D,2));
  
for i=1:budget
  disp(sprintf(' coverageNemhauserWolsey: finding element %d %s',i,datestr(now)));
  j = getNextIndexCoverageNemhauserWolsey(D,inds,covered);
  ind = [ind;inds(j)];
  covered = ((covered+D(inds(j),:))>0);
  inds = [inds(1:j-1) inds(j+1:end)];
end

end

function j = getNextIndexCoverageNemhauserWolsey(D,inds,covered)
  bestVal = -Inf;
  j = -1;
  for i=1:length(inds)
    val = sum((covered+D(inds(i),:))>0);
    if(val>bestVal)
      j = i;
      bestVal = val;
    end
  end
end
