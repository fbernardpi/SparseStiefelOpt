function [P] = v2p(v,c)
% V2P converts a vector representing a permutation into a permutation matrix
%       The function handles also partial permutations
% ----- Input:
%       v: vector representing a permutation
%       c: permutation size
% ----- Output:
%       P: (partial) permutation matrix
% ----- Author:
%       Andrea Fusiello 

n = length(v);
if nargin < 2
    c = n;
end
i = (1:n);
i(v==0) = [];
v(v==0) = [];
P = sparse(i,v,1,n,c);

% double stochastic
assert(all(sum(P,1)<=1) && all(sum(P,2)<=1))

end

