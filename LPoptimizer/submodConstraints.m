% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918


function A = submodConstraints(p)
    % A*x <= 0 is equivalent to x being a vector encoding a submodular
    % function

    % from the standard definitions, this has the least number of
    % constraints when all constraints are added naively:
    % A function is submodular iff for every B subset V and x1,x2 in V \ B
    % f(B U {x1}) + f(B U {x2}) >= f(B U {x1,x2}) + f(B)
    % -1*f(B U {x1}) + (-1)*f(B U {x2}) + 1*f(B U {x1,x2}) + 1*f(B) <= 0
  if(nargin<1)
    p = 5;
  end

  A = [];
    for i=1:2^p
        B = itovec(i-1,p);
        ind = find(B==0);
        for j=1:length(ind)-1
            for k = j+1:length(ind)
                v = zeros(1,2^p);
                v(i) = 1; % 1*f(B)
                Btilde = B;
                Btilde(ind(j)) = 1;
                v(vectoi(Btilde)+1) = -1; % -1*f(B U {x1})
                Btilde(ind(k)) = 1;
                v(vectoi(Btilde)+1) = 1; % 1*f(B U {x1,x2})
                Btilde = B;
                Btilde(ind(k)) = 1;
                v(vectoi(Btilde)+1) = -1; % -1*f(B U {x2})
                A = [A;v];
            end
        end
    end

end

