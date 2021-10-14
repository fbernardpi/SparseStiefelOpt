% 
% This function is a wrapper to apply the method in [1] to permutation 
% synchronisation. If you use this in your work, you are required to cite [1].
%
% [1] F. Bernard, D. Cremers, J. Thunberg. Sparse Quadratic Optimisation over 
% the Stiefel Manifold with Application to Permutation Synchronisation.
% NeurIPS 2021
%
% Input:   W is a matrix of size m x m, which comprises of k x k block
%              matrices that represent pairwise matchings, where the 
%              dimension of the block at position (i,j) is 
%              dimVector(i) x dimVector(j)
%          dimVector denotes the size of the blocks of W, where
%              sum(dimVector) = m
%          d is the number of columns of the output U (also called universe
%              size in permutation synchronisation problems)
%          vis (0 or 1) to enable visualisation
% Output:  Wout is a matrix of size m x m and contains cycle-consistent pairwise matchings
%          Uproj is a matrix of size m x d, so that Wout = Uproj*Uproj'
%
%
% Author and copyright: Florian Bernard (f.bernardpi@gmail.com)
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [Wout, Uproj] = SparseStiefelSync(W, dimVector, d, vis)

    U = SparseStiefelOpt(W, d, vis);
    Uproj = projectOntoPartialPermBlockwise(U, dimVector, [], 0);
    Wout = Uproj*Uproj';
  
end

