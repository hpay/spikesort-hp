% LOAD IN SORTED EPHYS DIRECTORIES FROM BIRD MASTER FOLDER
function [gm, idx,labels] = GMM_make(all_dir, option_only_good, n_clusters, label_by)
% label_by: how clusters are named E vs I (default 'rate')
%   'rate'  = lower mean log rate is E (original)
%   'width' = wider mean trough-to-peak duration is E (narrow-spiking = I)
% idx:
%   1 = ind_exc
%   2 = ind_inh
%   3 = ind_noise
% See this link on outputs: https://www.mathworks.com/help/stats/gmdistribution.cluster.html

%% Get waveforms from each session
[stats, goodIDs, mxWF, tWF] = GMM_stats(all_dir, option_only_good);

% Train only on negative-trough units: trough features don't apply when the peak is bigger
% (GMM_apply labels those units unknown)
pos_peak = max(mxWF,[],2) > abs(min(mxWF,[],2));
fprintf('GMM_make: excluding %d of %d units with peak > trough from training\n', nnz(pos_peak), numel(pos_peak))
stats = stats(~pos_peak,:);

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
if ~exist('label_by','var') || isempty(label_by); label_by = 'rate'; end
switch label_by
    case 'rate';  [~, order] = sort(mus.rate_log10);           % E first = lowest rate
    case 'width'; [~, order] = sort(mus.duration, 'descend');  % E first = widest spike
    otherwise; error('label_by must be ''rate'' or ''width''')
end
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

