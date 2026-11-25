%% Run the openEphys pipeline (runSingleKilosortOpenEphys) on the sessions listed in DataNotebook-MH.xlsx
% Sessions: the Folder column of the notebook's Ephys sheet. By default only sessions already sorted locally
% (the sort folder has spike_times.npy); with sort_new = true the others are also copied from the server and
% sorted (~50 GB and a Kilosort 4 sort each). Pipeline settings (cell type model, step 7 cutoffs, ...) are
% runSingleKilosortOpenEphys's defaults unless given in pipeline_opts.
% Run from the scripts folder (runSingleKilosortOpenEphys is here).

notebook = 'D:\dataMH\DataNotebook-MH.xlsx';
data_root = 'D:\dataMH';
ks_subdir = 'kilosort4_mirrored_interleaved_nodrift';
sort_new = false;            % true: also copy and sort sessions that aren't sorted locally yet
pipeline_opts = {};          % name-value settings for runSingleKilosortOpenEphys, e.g. {'min_type_p', 0.8, 'save_png', true}

T = readtable(notebook, 'Sheet', 'Ephys', 'VariableNamingRule', 'preserve');
sessions = T.Folder(~cellfun(@isempty, T.Folder));
is_sorted = cellfun(@(s) isfile(fullfile(data_root, s, ks_subdir, 'spike_times.npy')), sessions);
if ~sort_new
    fprintf('Not sorted locally, skipped: %s\n', strjoin(sessions(~is_sorted), ', '))
    sessions = sessions(is_sorted);
end
fprintf('Running %d sessions: %s\n', numel(sessions), strjoin(sessions, ', '))
for i = 1:numel(sessions)
    fprintf('\n===== %s =====\n', sessions{i})
    runSingleKilosortOpenEphys(sessions{i}, ks_subdir, 'data_root', data_root, pipeline_opts{:});
end
