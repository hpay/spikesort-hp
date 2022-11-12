% LOAD IN SORTED EPHYS DIRECTORIES FROM BIRD MASTER FOLDER
function [gm, idx,nlogL,P,logpdf,d2] = makeGMM_celltypes(all_dir, option_only_good, n_clusters)
% See this link on outputs: https://www.mathworks.com/help/stats/gmdistribution.cluster.html

%% Get waveforms from each session
[stats, goodIDs, mxWF, tWF] = getGMM_stats(all_dir, option_only_good);

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
[idx,nlogL,P,logpdf] = gm.cluster(X);

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

% Re-run GMM with new initialization: 1 = E, 2 = I, 3 = unknown
% My previous notation: % 0 = I, 1 = E, NaN = unknown
cov_type = 'full'; % SC uses diagonal but I find clusters are not always round
gm = fitgmdist(X,n_clusters,'CovarianceType',cov_type,'Start',ind_start);
[idx,nlogL,P,logpdf,d2] = gm.cluster(X);
P_plot = log(P(sub2ind(size(P),(1:length(idx))', idx)));
d2_plot = log(d2(sub2ind(size(d2),(1:length(idx))', idx)));
if n_clusters==2
    max_d2 = 2.5;
    idx(d2_plot>max_d2) = 3;
end

% Plot results
fs = 1/mean(diff(tWF));
figure; subplot(2,3,1);
scatter3( stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10, 15,stats.pk_deriv, 'filled');  axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
c=colorbar('north'); c.Label.String = 'Ratio peak:trough deriv.';

subplot(2,3,2)
scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,idx,'filled'); axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
subplot(2,3,3)
scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,d2_plot,'filled'); axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
c=colorbar('north'); 

%% Plot waveforms
subplot(2,3,4); plot(tWF,zscore(mxWF(idx==1, :)'),'k'); axis square; title('Excitatory')
hold on,plot(tWF,mean(zscore(mxWF(idx==1,:)'),2))
subplot(2,3,5); plot(tWF,zscore(mxWF(idx==2, :)'),'k'); axis square; title('Inhibitory')
if n_clusters==3
hold on,plot(tWF,mean(zscore(mxWF(idx==2,:)'),2))
subplot(2,3,6); plot(tWF,zscore(mxWF(idx==3, :)'),'k'); axis square; title('unknown')
end
