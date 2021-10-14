%
% Author and copyright: Florian Bernard (f.bernardpi@gmail.com)
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
function Z = computeZviaLineSearch(U, p)
    if ( ~exist('p', 'var') )
        p = 3;
    end
    d = size(U,2);
    
    
    try
        manifold = stiefelfactory(d,d,1);
    catch
        %% Setup ManOpt, which is used for linesearch
        addpath(genpath('./manopt'));
        
        manifold = stiefelfactory(d,d,1);
    end
    
    problem.M = manifold;

  
    % maximise the cost function trace((evec*Q)'*((evec*Q).^2)) s.t. Q'Q = I
    % (ManOpt minimises, so we use minus)
    problem.cost  = @(Q) -trace((U*Q)'*((U*Q).^(p-1)));
    problem.egrad = @(Q) -p*U'*((U*Q).^(p-1));
    %     checkgradient(problem);

    X = eye(d);

    [cost, grad] = getCostGrad(problem, X);
    % Pick the descent direction as minus the gradient
    desc_dir = problem.M.lincomb(X, -1, grad);
    gradnorm = problem.M.norm(X, grad);

    % Execute the line search
    [~, Z] = linesearch(problem, X, desc_dir, cost, -gradnorm^2);
    
end