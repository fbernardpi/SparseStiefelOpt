% Load data. The data in the file house_sync_data is created from the
% stereo matching code of Pachauri et al., NIPS, 2013 
% ( http://pages.cs.wisc.edu/~pachauri/perm-sync/ )
%
% It contains the following variables:
%      Wnoise          The (noisy) matrix of pairwise permutations
%      m               The number of points observed in each image
%      k               The total number of images
%      data            The keypoint coordinates in each image
%
% The ground truth is given by the identity pairwise matching matrix, i.e.
% Pgt = kron(ones(k),speye(m));
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
load house_sync_data

addpath(genpath(pwd));


%%
dimVector = repmat(m,k,1);
Wcell = mat2cell(Wnoise, dimVector, dimVector);


allErrors = [];
allObjectives = [];
allRuntimes = [];
allFscores = [];
allMethods = {'MatchEig', 'MatchALS', 'Spectral', 'NmfSync', 'Ours'};

% If you want to run MatchALS it is recommended to install MOSEK (also see
% README.txt)
% addpath(genpath('~/code/git/code/src/toolboxes/mosek/9.1')) 
runMatchALS = 0; % disabled by default for runtime reasons

% number of objects that are used to create gradually increasing permutation 
% synchronisation problems
kValues = 20:111; 
for i=1:numel(kValues)
    % note that this variable only stores the obtained matchings from the
    % most recent iteration (due to memory reasons); it is later used to
    % visualse these results (see further below)
    allP = {}; 
    
    % number of objects we synchronise
    currK = kValues(i);
    
    disp(['Processing k = ' num2str(currK)]);
    
    % subsample the original pairwise matchings to create a sequence of
    % permutation synchronisation problems with increasing size
    subIdx = round(linspace(1,k,currK));
    
    % noisy pairwise matching matrix
    currWnoise = cell2mat(Wcell(subIdx,subIdx));
    
    % vector containing the number of points for each image
    currDimVector = dimVector(subIdx);
    
    % ground truth
    currPgt = kron(ones(currK),speye(m));
    
    % synchronisation objective that we want to maximise
    obj = @(U) trace(U'*currWnoise*U);
    
    %% MatchEig
    s = tic();
    P_matchEig = MatchEIG(currWnoise,m,currK,currDimVector, 0);
    timeMatchEig = toc(s);
    objMatchEig = nan; % not cycle-consistent, so there does not exist a binary U matrix that would allow us to call obj()
    fscoreMatchEig = fscore(currPgt, P_matchEig);
    allP{end+1} = P_matchEig;
    
    %% MatchALS
    if ( runMatchALS )
        s = tic();
        Z_matchALS = mmatch_CVX_ALS(currWnoise,currDimVector,'pselect',0.7, ...
            'maxiter',100,'univsize',2*m,'verbose',false);
        timeMatchAls = toc(s);
        fscoreMatchAls = fscore(currPgt, Z_matchALS);
        allP{end+1} = Z_matchALS;
    else
        timeMatchAls = nan;
        fscoreMatchAls = nan; %fscore(currPgt, Z_matchALS);
        allP{end+1} = [];
    end
    objMatchAls = nan; % not cycle-consistent, so there does not exist a binary U matrix that would allow us to call obj()
    
    
    %% Spectral
    s = tic();
    [Pspectral,~,~,U_spectral] = mmatch_spectral(currWnoise,currDimVector,m); 
    
    timeMatchSpectral = toc(s);
    objSpectral =  obj(U_spectral);
    fscoreSpectral = fscore(currPgt, Pspectral);
    allP{end+1} = Pspectral;
    
    %% NmfSync
    s = tic();
    [ZnmfSync, UnmfSync] = nmfSync(full(currWnoise), currDimVector, m, []);
    timeMatchNmfSync = toc(s);    
    objNmfSync = obj(UnmfSync);
    fscoreNmfSync = fscore(currPgt, ZnmfSync);
    allP{end+1} = ZnmfSync;
    
    %% Ours
    s = tic();
    [P_our, U_our] = SparseStiefelSync(currWnoise, currDimVector, m, 0);
    timeOur = toc(s);
    objOur = obj(U_our);
    fscoreOur = fscore(currPgt, P_our);
    allP{end+1} = P_our;
    
    %% store all results
    allObjectives(end+1,:) = [objMatchEig objMatchAls objSpectral objNmfSync objOur]./currK^2; 
    allRuntimes(end+1,:) = [timeMatchEig timeMatchAls timeMatchSpectral timeMatchNmfSync timeOur];
    allFscores(end+1,:) = [fscoreMatchEig fscoreMatchAls fscoreSpectral fscoreNmfSync fscoreOur];
end


%% create plots
lineWidth = 2;

figure('color', 'w', 'position', [100 100 1200 300]);

linetypes = {'-.','-.','-', '-','-'};


nMethods = numel(allMethods);
colorMap = [0.650980392156863,0.807843137254902,0.890196078431373;0.121568627450980,0.470588235294118,0.705882352941177;0.698039215686275,0.874509803921569,0.541176470588235;0.200000000000000,0.627450980392157,0.172549019607843;0.984313725490196,0.603921568627451,0.600000000000000];

% fscore
subplot 131
hold on;
for jj=1:nMethods
    plot(kValues,allFscores(:,jj)', ...
        'linewidth', lineWidth, ...
        'color', colorMap(jj,:), 'linestyle', linetypes{jj});

end
legend(allMethods, 'orientation', 'horizontal');
xlim([kValues(1) kValues(end)]);
xlabel('k');
ylabel('fscore');


% objectives
subplot 132
hold on;
for jj=3:nMethods
    plot(kValues,allObjectives(:,jj)', ...
        'linewidth', lineWidth, ...
        'color', colorMap(jj,:), 'linestyle', linetypes{jj});

end

xlim([kValues(1) kValues(end)]);
xlabel('k');
ylabel('objective');


% runtime
subplot 133
hold on;
for jj=1:nMethods
    plot(kValues,allRuntimes(:,jj)', ...
        'linewidth', lineWidth, ...
        'color', colorMap(jj,:), 'linestyle', linetypes{jj});

end

xlim([kValues(1) kValues(end)]);
xlabel('k');
ylabel('runtime [s]');


    

%% show result images
imFolder = 'house_images/';

lw = 2;
sz = 200;
  
% add new rows for any pair of image indices (between 1 and 111) to visualise
subIdx = [1 111];
for i=1:size(subIdx,1)
    %%
    I1 = imread([imFolder 'house.seq' num2str(subIdx(i,1)-1) '.png']);
    I2 = imread([imFolder 'house.seq' num2str(subIdx(i,2)-1) '.png']);
    [r,c] = size(I1);
    
    pts1 = data{subIdx(i,1)};
    pts2 = data{subIdx(i,2)};
    
    % crop images to allow for a tighter visualisation
    I1LeftOffset = 140;
    I1RightOffset = 510;
    I2LeftOffset = 80;
    I2RightOffset = 360;
    
    I2(:,1:I2LeftOffset) = [];
    I2(:,I2RightOffset:end) = [];
    I1(:,I1RightOffset:end) = [];
    
    pts1(:,1) = pts1(:,1) - I1LeftOffset;
    pts2(:,1) = pts2(:,1) - I1LeftOffset + I1RightOffset - I2LeftOffset;

    
    II = [I1, I2];
    II(:,1:I1LeftOffset) = [];

    cols = lines(size(pts1,1));
    
    for methodIdx=1:numel(allMethods)
        if ( ~isempty(allP{methodIdx}) )
            matchingGt = 1:size(pts1,1);
            
            iidx = (1+(subIdx(i,1)-1)*m):subIdx(i,1)*m;
            jidx = (1+(subIdx(i,2)-1)*m):subIdx(i,2)*m;
            
            currP = allP{methodIdx}(iidx,jidx);
            currMatching = zeros(1,m);
            [i1,i2] = find(currP);
            currMatching(i2) = i1;
            
            correctIdx = find(matchingGt==currMatching);
            wrongIdx = find(matchingGt~=currMatching & currMatching~=0);
            missingIdx = find(currMatching==0);
            
            
            figure
            
            imshow(II)
            hold on
            
            scatter(pts1(:,1), pts1(:,2),sz, cols, 'filled', 'markeredgecolor', 'k', 'linewidth', 2)
            scatter(pts2(:,1), pts2(:,2),sz, cols, 'filled', 'markeredgecolor', 'k', 'linewidth', 2)
            
            line([pts1(correctIdx,1), pts2(currMatching(correctIdx),1)]', [pts1((correctIdx),2), pts2(currMatching(correctIdx),2)]', 'color', 'g', 'linewidth', lw);
            line([pts1((wrongIdx),1), pts2(currMatching(wrongIdx),1)]', [pts1((wrongIdx),2), pts2(currMatching(wrongIdx),2)]', 'color', 'r', 'linewidth', lw);
            
            drawnow;
            
            title(allMethods{methodIdx})
        end
    end
end




     