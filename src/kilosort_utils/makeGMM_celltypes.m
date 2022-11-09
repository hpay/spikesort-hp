% LOAD IN SORTED EPHYS DIRECTORIES FROM BIRD MASTER FOLDER
function gm = makeGMM_celltypes(allDir, option_only_good)

%% Get waveforms from each session
for ii = 1:length(allDir)
    ks_dir = fullfile(allDir{ii},'kilosort2_output');
    wvStruct = getfield(load(fullfile(ks_dir,'waveformStruct.mat')),'wvStruct');
    
    % Only keep "good" clusters. If sorted with Phy, this will be the final sorting
    if ~exist('option_only_good','var') || option_only_good
        [unit_ID, cluster_labels] = get_phy_cluster_labels(ks_dir); % 0 = noise, 1 = 'mua', 2 = 'good'
        mask = (cluster_labels==2);
        wvStruct.mxWF = wvStruct.mxWF(mask,:);
        wvStruct.max_site = wvStruct.max_site(mask,:);
        wvStruct.pcWF = wvStruct.pcWF(mask,:);
        wvStruct.meanRate = wvStruct.meanRate(mask,:);
        wvStruct.waveFormsMean = wvStruct.waveFormsMean(:,:,mask);
        wvStruct.goodIDs = wvStruct.goodIDs(mask);
        wvStruct.goodLabels = wvStruct.goodLabels(mask);
        wvStruct.medISI = wvStruct.medISI(mask);
        wvStruct.contam = wvStruct.contam(mask);
    end
    
    % Add mxWF_d
    wvStruct.mxWF_d = diff(wvStruct.mxWF')';
    theseWaveforms(ii) = wvStruct;
end


%% Compute stats
meanRate = cat(1,theseWaveforms.meanRate);
mxWF = cat(1,theseWaveforms.mxWF);
mxWF_d = cat(1,theseWaveforms.mxWF_d);

numUnits = length(meanRate);
stats = array2table(NaN(numUnits,4)); stats.Properties.VariableNames = {'rate_log10','duration','asymmetry','pk_deriv'};
for nUnit = 1:numUnits
    % get trough
    [tr, trInd] = min(mxWF(nUnit,:));
    
    % get peak (after trough)
    [pk, pkInd] = max(mxWF(nUnit, trInd:end));
    pkInd = pkInd + trInd - 1;
    
    % get pre-peak (before trough)
    pk_pre = max(mxWF(nUnit, 1:trInd));
    
    % get derivative trough and peak
    [tr_d,trInd_d] = min(mxWF_d(nUnit,:));
    [pk_d,pkInd_d] = max(mxWF_d(nUnit,:));
    
    % concatenate the following statS  
    stats.rate_log10(nUnit) = log10(meanRate(nUnit)); 
    stats.duration(nUnit) =  pkInd-trInd;            % duration
    stats.asymmetry(nUnit) = (pk-pk_pre)./(pk+pk_pre);      % asymmetry
    stats.pk_deriv(nUnit) = log10(abs(pk_d./tr_d));  % pk derivative ratio
end

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
fs = size(wvStruct.mxWF,2)/wvStruct.spkDur;

figure; subplot(1,2,1);
scatter3( stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10, 15,stats.pk_deriv, 'filled');  axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
c=colorbar('north'); c.Label.String = 'Ratio peak:trough deriv.';

subplot(1,2,2)
scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,idx,'filled'); axis square; axis vis3d
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')

%% Plot waveforms
x = linspace(-theseWaveforms(1).spkOffset,-theseWaveforms(1).spkOffset+theseWaveforms(1).spkDur,...
    fs*theseWaveforms(1).spkDur);
figure;
subplot(1,3,1); plot(x,zscore(mxWF(idx==1, :)'),'k'); axis square; title('Excitatory')
hold on,plot(x,mean(zscore(mxWF(idx==1,:)'),2))
subplot(1,3,2); plot(x,zscore(mxWF(idx==2, :)'),'k'); axis square; title('Inhibitory')
hold on,plot(x,mean(zscore(mxWF(idx==2,:)'),2))
subplot(1,3,3); plot(x,zscore(mxWF(idx==3, :)'),'k'); axis square; title('unknown')

