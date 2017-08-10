

function A = setfnMaxRandomGreedy(g,p,init1,numRandomSeeds);
% matthew.blaschko@esat.kuleuven.be
% copyright 2016
%
% finds a local maximum of an arbitrary set function g by steepest ascent
% from a set of random seeds.  It also will always try searching from the
% empty set
%
% inputs:
%  g - a set function
%  p - the size of the base set, |V|
%  init1 - an imputed point to maximize from
%  numRandomSeeds - number of random binary vectors of length p to start
%                   from, by default 10

if(nargin<4)
    numRandomSeeds = 10;
end

init = zeros(p,1);

A = setfnMaxGreedy(g,init,p);

if(nargin>=3)
  currA = setfnMaxGreedy(g,init1,p);
  if(g(find(currA))>g(find(A)))
    A = currA;
  end
end

for i=1:numRandomSeeds
    init = round(rand(p,1));
    currA = setfnMaxGreedy(g,init,p);
    if(g(find(currA))>g(find(A)))
        A = currA;
    end
end

A = find(A);

end

function A = setfnMaxGreedy(g,A,p)

origA = A;

bestScore = g(find(A));

for i=1:p
    currA = origA;
    currA(i) = 1 - currA(i);
    currScore = g(find(currA));
    if(currScore>bestScore)
        bestScore = currScore;
        A = currA;
    end
end

if(norm(A-origA)>0)
    A = setfnMaxGreedy(g,A,p);
end

end

