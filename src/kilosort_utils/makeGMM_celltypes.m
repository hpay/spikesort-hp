% LOAD IN SORTED EPHYS DIRECTORIES FROM BIRD MASTER FOLDER
function gm = makeGMM_celltypes(all_dir, option_only_good)

%% Get waveforms from each session
[stats, goodIDs, mxWF, tWF] = getGMM_stats(all_dir, option_only_good);

% Exclude very low rate cells from the template creation
meanRateThresh = 0.05;  % min rate threshold
mask = stats.rate_log10 < log10(meanRateThresh); 
stats(mask,:) = [];
mxWF(mask,:) = [];

%% Run GMM
rng(2)
X = table2array(stats);
gm = fitgmdist(X,3,'Replicates',25,'CovarianceType','diagonal');
[idx,nlogL,P,logpdf] = gm.cluster(X);

% Identify clusters
mus = array2table(gm.mu); mus.Properties.VariableNames = stats.Properties.VariableNames;
sigmas = array2table(squeeze(gm.Sigma)'); sigmas.Properties.VariableNames = stats.Properties.VariableNames;
sigmas_test = table2array(sigmas(:,{'rate_log10','duration'}));
[~, ind_noise] = max(mean(sigmas_test./mean(sigmas_test),2));
[~, order] = sort(mus.rate_log10);
order(order==ind_noise)=[];
ind_exc = order(1);
ind_inh  = find(~ismember(1:3, [ind_exc ind_noise]));
ind_start = zeros(height(stats),1);
ind_start(idx==ind_exc) = 1;
ind_start(idx==ind_inh) = 2;
ind_start(idx==ind_noise) = 3;

% Re-run GMM with new initialization: 1 = E, 2 = I, 3 = unknown
% My previous notation: % 0 = I, 1 = E, NaN = unknown
gm = fitgmdist(X,3,'CovarianceType','diagonal','Start',ind_start);
[idx,nlogL,P,logpdf] = gm.cluster(X);

% Plot results
fs = 1/mean(diff(tWF));
figure; subplot(1,2,1);
scatter3( stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10, 15,stats.pk_deriv, 'filled');  axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
c=colorbar('north'); c.Label.String = 'Ratio peak:trough deriv.';

subplot(1,2,2)
scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,idx,'filled'); axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')

%% Plot waveforms

figure;
subplot(1,3,1); plot(tWF,zscore(mxWF(idx==1, :)'),'k'); axis square; title('Excitatory')
hold on,plot(tWF,mean(zscore(mxWF(idx==1,:)'),2))
subplot(1,3,2); plot(tWF,zscore(mxWF(idx==2, :)'),'k'); axis square; title('Inhibitory')
hold on,plot(tWF,mean(zscore(mxWF(idx==2,:)'),2))
subplot(1,3,3); plot(tWF,zscore(mxWF(idx==3, :)'),'k'); axis square; title('unknown')

