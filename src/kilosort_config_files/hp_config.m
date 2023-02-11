function ops = hp_config(ops)

if ~exist('ops','var') || isempty(ops) || ~isfield(ops,'chanMap')
    ops.Nchan = 64;
    ops.NchanTOT = 64;
    ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file    
else
    temp = load(ops.chanMap);    
    ops.Nchan = numel(temp.chanMap);
    ops.NchanTOT = numel(temp.chanMap);    
end

% Temporary scratch space, ideally on SSD?
ops.fproc   = 'C:\temp\tmp_wh.dat';

% time range to sort
ops.trange    = [0 Inf]; 

% sample rate
ops.fs = 30000;  

% frequency for high pass filtering 
ops.fshigh = 300;           % (SC changed from 150 to 300)

% minimum firing rate on a "good" channel (used during pre-processing) (0 to skip)
ops.minfr_goodchannels = 0; % (SC changed from 1/50 to 0). 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 10;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50;  % (SC changed from 1/50 to 1/100)

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 

% SC: subtract channel mean before CAR, should not be needed for our probes
ops.chMeanSub = 0; 

% SC: Use CAR (actually median)
ops.CAR = 1; 

% artifact detection
ops.detectArtifacts = 0; % (1) artifact detection, uses parameters below to eliminate single-sample positive going blips
ops.artifactStdThresh = [40, 8];% ([40,8]) in (robust) std of derivative, for global and single-channel respectively
ops.artifactRatioThresh = 1/3; 	% (1/3) abs. value of sum of thresh and thresh+1 must be below this (for single-channel)
ops.blankingSmWin = 11;  		% smooth voltage derivates before peak finding
ops.blankingHalfWidth = round(0.1 * ops.fs); % window around peaks to blank data

% LFP
ops.makeLFP = 0;
ops.dsFac = 32; % this must cleanly divide into ops.NT without remainder
ops.fsLFP = ops.fs/ops.dsFac; %1,000hz
ops.lpFreq = 300; % 300hz lowpass cutoff

%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;       % (-6) spike threshold in standard deviations 
ops.reorder         = 0;        % (SC changed from 1 to 0, but in run_batch_kilosort script) whether to reorder batches for drift correction. % *SC/HP: set to 0 to order batches linearly by time
ops.nskip           = 20;       % (SC changed from 25) how many batches to skip for determining spike PCs *HP: why did SC change?***

ops.GPU                 = 1;    % has to be 1, no CPU version yet
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 12;    % (SC changed from 4) % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 128;  % (SC changed from 64) % samples of symmetrical buffer for whitening and spike detection
ops.NT = 8*64*1024 + ops.ntbuff;% (SC changed from 64*1024+ ops.ntbuff) % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 64*1024/30k = 2.1845 s, 8*64*1024/30k = 17.5 s
ops.whiteningRange      = 32;   % number of channels to use for whitening each channel
ops.nSkipCov            = 20;   % (SC changed from 25) % compute whitening matrix from every N-th batch *HP: why?
ops.scaleproc           = 200;  % int16 scaling of whitened data 
ops.nPCs                = 3;    % how many PCs to project the spikes into
ops.useRAM              = 0;    % not yet available

ops.nt0 	            = 61;   % number of time samples for the templates (has to be <=81 due to GPU shared memory)
ops.nt0min = floor(ops.nt0/3);  % time sample where the negative peak should be aligned
