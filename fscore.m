%
% Function to compute the fscore
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
function [fscore,precision,recall] = fscore(Wgt, Win)
    Wgt = triu(Wgt,1);
    Win = triu(Win,1);

    TP = nnz(Wgt & Win);
    FP = nnz(~Wgt & Win);
    FN = nnz(Wgt & ~Win);

    precision = TP./(TP+FP);
    recall = TP./(TP+FN);

    fscore = (2*precision*recall)./(precision + recall);
end