%% Function to do sorting that takes in root data directory and options
function run_single_kilosort(rootDir, saveDir, ops)

% find the binary file
fprintf('Looking for data inside %s \n', rootDir)
ops.fbinary = fullfile(rootDir, 'amplifier.dat'); % the raw data binary file is in this folder

% run main algo
rez = preprocessDataSub(ops);    % preprocess data to create temp_wh.dat
rez = clusterSingleBatches(rez); % time-reordering as a function of drift
rez = learnAndSolve8b(rez);      % main tracking and template matching algorithm

% OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort2/issues/29
rez = remove_ks2_duplicate_spikes(rez,'overlap_s',5e-5,'channel_separation_um',30);

% Final processing
rez = find_merges(rez, 1);      % final merges
rez = splitAllClusters(rez, 1); % final splits by SVD
rez = splitAllClusters(rez, 0); % final splits by amplitudes
rez = set_cutoff(rez); % decide on cutoff
fprintf('found %d good units \n', sum(rez.good>0))

% Write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, saveDir);


%% If you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
fprintf('Saving results to Rez  \n')
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez
fname = fullfile(saveDir, 'rez.mat');
save(fname, 'rez', '-v7.3');
end