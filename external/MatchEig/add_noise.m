function Z = add_noise(Z,n_err)
% ADD_NOISE introduces mismatches in the permutation matrix Z
% ----- Input:
%       Z: ground truth permutation matrix
%       n_err: number of mismatches 
% ----- Output:
%       Z: noisy permutation matrix
% ----- Authors:
%       Eleonora Maset, Federica Arrigoni and Andrea Fusiello, 2017  

d = size(Z,1);
c = size(Z,2);
count_err = 0;  % # mismatches introduced

while count_err < n_err
    
    matchZ = find(sum(Z,2) == 1);	
    no_match_row = find(sum(Z,2) == 0);
    no_match_col = find(sum(Z,1) == 0);
    
    p = rand;   % select the type of mismatch
    
    % no common features between the two images
    if isempty(matchZ)
        p = 0.9;
    end
    
    % identify a total permutation
    if sum(sum(Z)) == min(d,c)
        p = 1/3;
    end
    
    if p < 1/3 && n_err-count_err > 1
        % switch two matches
        idx_match_1 = matchZ(randi(length(matchZ)));
        match_1 = Z(idx_match_1,:);
        idx_match_2 = randi(d);
        while idx_match_2 == idx_match_1
            idx_match_2 = randi(d);
        end
        match_2 = Z(idx_match_2,:);
        Z(idx_match_1,:) = match_2;
        Z(idx_match_2,:) = match_1;
        
        count_err = count_err + 2;
        
    elseif p >= 1/3 && p < 0.9
        % remove a true match
        idx_match = matchZ(randi(length(matchZ)));
        Z(idx_match,:) = zeros(1,size(Z,2));
        
        count_err = count_err+1;
        
    else
        % add a false match
        idx_match_1 = no_match_row(randi(length(no_match_row)));
        idx_match_2 = no_match_col(randi(length(no_match_col)));
        Z(idx_match_1,idx_match_2) = 1;
        
        count_err = count_err+1;
        
    end   
end
end