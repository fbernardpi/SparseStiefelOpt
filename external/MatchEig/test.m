% Demonstration of the multi-view matching algorithm MatchEIG
% Paper: Maset E., Arrigoni F., Fusiello A.: Practical and Efficient
%        Multi-view matching - ICCV'17

clear all

% problem setting
n = 30;             % # views
d = 100;            % # features
obs_ratio = 0.8;    % observation ratio 
err_rate = 0.2;     % input error (ratio of observations corrupted)

% parameter setting (MatchEIG) 
thresh = 0.25;

% Create synthetic data
P = [];
dimPerm = zeros(n,1);

for i = 1:n
    v = zeros(d,1);
    while sum(v) == 0
        v = randperm(d);
        v(rand(d,1)>obs_ratio) = 0;
    end
    p = v2p(v);
    iz = find(sum(p,2) == 0);
    p(iz,:) = [];
    dimPerm(i) = d-length(iz);  % # features in image i
    P = [P; p]; % p is the ground truth absolute permutation   
end

% Compute ground truth relative permutations
Z = P*P';
Z_gt = sparse(Z);

% Add noise to each relative permutation
cumDim = [0; cumsum(dimPerm(1:end-1))];
for i = 1:n
    for j = (i+1):n
        Zij = Z(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j));
        nf = min(size(Zij));
        n_err = round(nf*err_rate);   % # mismatches to introduce
        Znoise = add_noise(Zij,n_err);  
        % replace ground truth in Z with noisy relative permutations
        Z(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)) = Znoise;
        Z(1+cumDim(j):cumDim(j)+dimPerm(j),1+cumDim(i):cumDim(i)+dimPerm(i)) = Znoise';			
    end
end
Z = sparse(Z);  % input noisy matrix

% Run Multi-view matching (MatchEIG)
tic
[Z_mvm] = MatchEIG(Z,d,n,dimPerm,thresh);
time_match = toc;

% Error evaluation (F-score)
TP = 0;     % true positive
FPN = 0;    % false positive + false negative
for i = 1:n
    for j = 1:n
        TP = TP + sum(sum(Z_mvm(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)) & Z_gt(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)))) ;
        FPN = FPN + sum(sum(Z_mvm(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)) | Z_gt(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)))) - sum(sum(Z_mvm(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j)) & Z_gt(1+cumDim(i):cumDim(i)+dimPerm(i),1+cumDim(j):cumDim(j)+dimPerm(j))));
    end
end
F_score = 2*TP/(2*TP+FPN);

% Print results
fprintf('Results: Fscore = %.2f \nTime = %.2f s\n',F_score,time_match)


