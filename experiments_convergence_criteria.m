%% A simple toy example that demonstrates two different convergence criteria
% (please note that ManOpt (https://www.manopt.org) needs to be installed to run this script)
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
clc

W = diag([1,1,1,.5,.5,.5]); 
d = 3;

% Note that any choice of p works, where p is an integer >= 3
p = 3;

U = SparseStiefelOpt2(W, d, 0, p, 'primary');
U

disp('This small toy problem terminates too early, since the convergence criterion ');
disp('above is purely based on the primary objective tr(U^T*W*U). If we additionally consider ');
disp('the secondary objective as convergence criterion (which in turn requires linesearch ');
disp('in order to guarantee the secondary objective is monotonically increasing), ');
disp('the result becomes sparse:');

U2 = SparseStiefelOpt2(W, d, 0, p, 'primaryAndSecondary');
U2

figure;
subplot 121
imagesc(U)
colorbar
title('tr(U^TWU) as sole convergence criterion');

subplot 122
imagesc(U2)
colorbar
title('tr(U^TWU) and sparsity as convergence criterion');