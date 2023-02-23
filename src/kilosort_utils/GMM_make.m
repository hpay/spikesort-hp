% LOAD IN SORTED EPHYS DIRECTORIES FROM BIRD MASTER FOLDER
function [gm, idx,labels] = GMM_make(all_dir, option_only_good, n_clusters)
% See this link on outputs: https://www.mathworks.com/help/stats/gmdistribution.cluster.html

%% Get waveforms from each session
[stats, goodIDs, mxWF, tWF] = GMM_stats(all_dir, option_only_good);

% Exclude very low rate cells from the template creation
% meanRateThresh = 0.05;  % min rate threshold
% mask = stats.rate_log10 < log10(meanRateThresh); 
% stats(mask,:) = [];
% mxWF(mask,:) = [];

if ~exist('n_clusters','var')
    n_clusters= 3;
end

%% Run GMM
rng(2)
X = table2array(stats);
gm = fitgmdist(X,n_clusters,'CovarianceType','diagonal','Replicates',25); % 
[idx,~,P,logpdf, d2] = gm.cluster(X);

% Identify clusters
mus = array2table(gm.mu); mus.Properties.VariableNames = stats.Properties.VariableNames;
[~, order] = sort(mus.rate_log10);
sigmas = array2table(squeeze(gm.Sigma)'); sigmas.Properties.VariableNames = stats.Properties.VariableNames;
sigmas_test = table2array(sigmas(:,{'rate_log10','duration'}));

if n_clusters==3
    [~, ind_noise] = max(mean(sigmas_test./mean(sigmas_test),2));
    order(order==ind_noise)=[];    
    ind_exc = order(1);
    ind_inh  = find(~ismember(1:3, [ind_exc ind_noise]));
    ind_start = zeros(height(stats),1);
    ind_start(idx==ind_exc) = 1;
    ind_start(idx==ind_inh) = 2;
    ind_start(idx==ind_noise) = 3;
else
    ind_exc = order(1);    
    ind_start = zeros(height(stats),1);
    ind_start(idx==ind_exc) = 1;
    ind_start(idx~=ind_exc) = 2;
end

%% Re-run GMM with new initialization: 1 = E, 2 = I, 3 = unknown
% My previous notation: % 0 = I, 1 = E, NaN = unknown
cov_type = 'full'; % SC uses diagonal but I find clusters are not always round
gm = fitgmdist(X,n_clusters,'CovarianceType',cov_type,'Start',ind_start);

plot_on = 1;
[idx,labels] = GMM_apply(gm, all_dir, option_only_good, plot_on);

end

