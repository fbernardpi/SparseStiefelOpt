function [X,V,info,Y] = mmatch_spectral(W,dimGroup,k, eigMode)

k = min(k,size(W,1));
t0 = tic;

if ( ~exist('eigMode', 'var') )
    eigMode = 'eigs';
end

switch eigMode
    case 'eigs'
        [V,~] = eigs(W,k,'la');
    case 'eig'
        [V,eval] = eig(full(W));
        V = V(:,end-k+1:end);
end

% rounding
if ( all(dimGroup == k) )
    % full matching case (as presented by Pachauri et al.)
    V1 = V(1:k,1:k);
    
    V = V*V1';
    Y = projectOntoPartialPermBlockwise(V, dimGroup);
else
    Y = rounding(V(:,1:k),dimGroup,0.5);
end
X = sparse(double(Y))*sparse(double(Y))';

info.time = toc(t0);
                        
                        
                

    
    


