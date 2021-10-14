% 
% This function implements the algorithm for sparse quadratic optimisation
% over the Stiefel manifold introduced in [1]. If you use this in your
% work, you are required to cite [1].
%
% [1] F. Bernard, D. Cremers, J. Thunberg. Sparse Quadratic Optimisation over 
% the Stiefel Manifold with Application to Permutation Synchronisation.
% NeurIPS 2021
%
% Input:   W is a matrix of size m x m
%          d is the number of columns of the output U
%          vis (0 or 1) to enable visualisation
% Output:  U is a sparse (see below) matrix of size m x d that is a global optimum to
%              max_U trace(U'*W*U) s.t. U'*U = eye(d)
%
% Sparsity is characterised in terms of the secondory objectve function 
% g(U) = U(:)'*(U(:).^2)
%
% Note that this is the function that was used in the main experiments of [1].
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
function U = SparseStiefelOpt(W, d, vis)
    maxIter = 200;
    epsilon = 1e-5;
    
    m = size(W,1);

    % make W symmetric (without changing the optimiser)
    W = .5*(W + W'); 

    % make W PSD (without changing the optimiser)
    lambdaMin = eigs(W, 1, 'sa');
    Wpsd = W - (lambdaMin-eps)*speye(m);

    % function to measure the binaryness
    binaryness = @(U) U(:)'*(U(:).^2);

    rng(123); % for reproducability
    U = randn(m,d);

    if ( ~exist('vis', 'var') )
        vis = 0;
    end

    
    obj = -inf;
    binarynessObj = -inf;
    for nIter=1:maxIter
        % compute the Z matrix
        H = U'*(U.^2);
        tH = H - H';
        tH = tH./norm(tH(:),inf);

        Z = eye(d) + tH;

        % matrix multiplication (as in the Orthogonal Iteration algorithm
        % but with Z additionally on the right)
        WUZ = Wpsd*(U*Z);

        % perform QR decomposition
        [U, R, ~] = qr(WUZ, 0);
        signDiag = spdiags(sign(diag(R)),0,d,d); % retraction QR
        U = U*signDiag;
        R = signDiag*R;

        if ( vis )
            figure(1);
            clf;

            subplot(1,3,1);
            imagesc(U);


            subplot(1,3,2);
            plot(diag(R));

            subplot(1,3,3);
            imagesc(R);
            colorbar
            
            drawnow;
        end

        obj(nIter+1) = trace(U'*W*U);
        binarynessObj(nIter+1) = binaryness(U);
        if ( vis )
            disp(['objective: ' num2str(obj(nIter+1), '%.2f') ...
                ', binaryness: ' num2str(binarynessObj(nIter+1), '%.2f')]);
        end

        if ( obj(nIter)/obj(nIter+1) >= 1 - epsilon  )
            disp(['Converged after ' num2str(nIter) ' iterations']);
            break;
        end
    end 
end

