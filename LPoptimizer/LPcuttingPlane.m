% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918

function [x,sparm,state] = LPcuttingPlane(sparm, oldstate)
    
    if (~isfield(sparm,'convergenceThreshold'))
        sparm.convergenceThreshold = 1e-6;
    end

    maxIterations = 500;

    if (exist('oldstate','var'))
        state = oldstate;
    else
        state = relaxedLP(); % initialize state
    end
    if(isfield(sparm,'budget'))
        state.budget = sparm.budget;
    end
    
    sparm.x = zeros(sparm.p,1);
    state.x = sparm.x;

    %for i=1:sparm.p
    %  A = i;
    %  sparm.m(i) = sparm.func(A);
    %end
    
    sparm.m = modularLowerBoundNull(sparm.func,sparm.p);

    
    minIterations = 2;
    numIterations = 0;

    bestPrimalObjective = Inf;
    
    while (((bestPrimalObjective - state.objective) > sparm.convergenceThreshold ...
            || minIterations>0) && numIterations < maxIterations)
        
        numIterations = numIterations + 1;
        minIterations = minIterations - 1;

        switch(sparm.formulationType)
            case 'margin'
                [phi, b] = computeMarginConstraint(sparm, sparm.func);
            case 'slack'
                [phi, b] = computeSlackConstraint(sparm, sparm.func);
            case 'slackmargin'
                [phi, b] = computeSlackMarginConstraint(sparm, sparm.func);
            otherwise
                error('forumulationType not set to a known case')
        end

        primalobjective = b - dot(state.x, phi);
        if (primalobjective < bestPrimalObjective)
            bestPrimalObjective = primalobjective;
            bestState = state;
            bestState.primalobjective = primalobjective;
        end
        

        fprintf([' %d primal objective: %f, best primal: %f, dual objective: %f, gap: %f\n'], ...
                   numIterations, primalobjective, bestPrimalObjective, state.objective, ...
                   (bestPrimalObjective - state.objective));

        state = relaxedLP(state, phi, b);
        sparm.x = state.x;
        model.x = state.x;
    end
    
    sparm.x = bestState.x;
    model.x = bestState.x;
    x = model.x;
    sparm.bestobjective = bestState.primalobjective;

end

function [phi, b] = computeSlackMarginConstraint(sparm,func)
    [phi, b, A] = computeMarginConstraint(sparm,func);
    [phi2,b2] = computeSlackConstraint(sparm,func,A);
    if(b2 - dot(sparm.x, phi2) > b - dot(sparm.x, phi))
        phi = phi2;
        b = b2;
    end
    [phi2,b2] = computeJegelkaConstraint(sparm,func,A);
    if(b2 - dot(sparm.x, phi2) > b - dot(sparm.x, phi))
        phi = phi2;
        b = b2;
    end
    [phi2,b2] = computeJ2constraint(sparm,func,A);
    if(b2 - dot(sparm.x, phi2) > b - dot(sparm.x, phi))
        phi = phi2;
        b = b2;
    end
end

function [phi, b, A] = computeMarginConstraint(sparm,func,A)
  if(~exist('A','var'))
   A = sfo_min_norm_point(@(A)(marginFunc(A,func,sparm.x,sparm.m)),[1:sparm.p]);
   %A = setfnMaxBruteForce(@(A)(-marginFunc(A,func,sparm.x,sparm.m)),sparm.p);
  end
%marginFunc(A,func,sparm.x,sparm.m)
%marginFunc(A2,func,sparm.x,sparm.m)
%assert(norm(marginFunc(A,func,sparm.x,sparm.m)-marginFunc(A2,func,sparm.x,sparm.m))<1e-9);
  tmp = zeros(sparm.p,1);
  tmp(A) = 1;
  phi = -(sparm.m' + tmp);
  b = func(A) - sum(sparm.m(A)) - length(A);
end

function val = marginFunc(A,func,x,m)
  A = unique(A);
  val = -(func(A) - sum(m(A)) - length(A) + sum(x(A)));
end

function [phi, b, A] = computeSlackConstraint(sparm,func,A)
  %isSubmodular(@(A)(slackFunc(A,func,sparm.x,sparm.m)),sparm.p) % not supermodular
  %A = setfnMaxRandomGreedy(@(A)(slackFunc(A,func,sparm.x,sparm.m)),sparm.p,round(sparm.x));
  if(~exist('A','var'))
      A = setfnMaxBruteForce(@(A)(slackFunc(A,func,sparm.x,sparm.m)),sparm.p);
  end
%A
  tmp = zeros(sparm.p,1);
  tmp(A) = 1;
  tmp = (func(A) - sum(sparm.m(A)))*tmp;
  phi = -(sparm.m' + tmp);
  b = (func(A) - sum(sparm.m(A))) *(1- length(A));
end

function val = slackFunc(A,func,x,m)
  val = (func(A) - sum(m(A)))*(1-length(A) + sum(x(A)));
end


