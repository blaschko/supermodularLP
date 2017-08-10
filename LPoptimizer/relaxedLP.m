% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918



%   Solves the problem
%
%    min_{x,xi} xi,  xi >= b_t - <a_t, x>  for t = 1, ..., T
%
%    0<=x<=1
%
%    subject to optional constraint sum x <= state.budget
%    usage:
%      state = relaxedLP();
%      state.budget = 15;
%      state = relaxedLP(state, a, b);
%      repeat until all desired constraints are added.

function state = relaxedLP(state, a, b)
    
    if nargin == 0
      state.a = [];
      state.b = [];
      state.Age = [];
      state.objective = -inf;
      state.x = [];
      return;
    end

    dimension = size(state.x,1);
    numNewConstraints = size(a,2);
    if isempty(state.a)
        state.a = zeros(dimension, 0);
        if(isfield(state,'budget'))
          state.harda = -ones(dimension,1); %add budget hard constraint \sum x <= state.budget
          state.hardb = -state.budget;
        else
          state.harda = [];
          state.hardb = [];
        end
    end;

    % add new constraints to the pool
    state.Age       = [state.Age ; 0];

    state.a = [state.a, a];
    state.b = [state.b, b];


    
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if(size(state.a,1)>0)
      if(~isOctave)
        [state.x,state.objective,exitflag,output,lambda] =...
	        linprog([1;zeros(size(state.a,1),1)],...
		              -[[zeros(size(state.harda,2),1);ones(size(state.a,2),1)] [state.harda';state.a']],-[state.hardb state.b],...
                  [],[],[-Inf;zeros(size(state.a,1),1)],...
                  [Inf;ones(size(state.a,1),1)]);
      else
        % CALL OCTAVE VERSION
        [state.x, state.objective, exitflag, EXTRA] =...
          glpk([1;zeros(size(state.a,1),1)],... % C
          -[[zeros(size(state.harda,2),1);ones(size(state.a,2),1)] [state.harda';state.a']],... %-[ones(size(state.a,2),1) state.a'],... % A
               -[state.hardb state.b],... %-state.b,... % B
               [-Inf;zeros(size(state.a,1),1)],... % LB
               [Inf;ones(size(state.a,1),1)],... % UP
               repmat('U',[1,length(state.b)+length(state.hardb)]),... % CTYPE
               repmat('C',[1,1+size(state.a,1)])... % VARTYPE
               );
        lambda.ineqlin = abs(EXTRA.lambda);
      end
    else
      state.x = [];
      state.objective = -Inf;
      lambda = [];
    end

    % remove idle variables
    state.Age = state.Age + 1;
    active = lambda.ineqlin > 1e-5;
    active = active(size(state.harda,2)+1:end);
    state.Age(active) = 0;
    keep = state.Age < Inf;

    state.Age = state.Age(keep);
    state.a = state.a(:,keep);
    state.b = state.b(keep);

    % update the model
    state.x = state.x(2:end);
