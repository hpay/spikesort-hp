%% Script for batch processing kilosort spike sorting
% and backing up recently collected data to server!
%
% Instructions: Update RECORDING_DEPTH_CHICK.xlsx then run this script from
% within the scripts folder
%
% Necessary column names in spreadsheet:
%   filename: name of the session directory
%   bad_chan: list any bad channels in this format: [1, 6, 10]. If none: []
%   probe: name of probe channel map {H5, H6}
%   manually_sorted: 1 if manually curated in phy, 0 otherwise
% 
% Currently set up to analyze data on local hard drive

Tall = readtable('C:\Users\Hannah\Dropbox\alab\Code\project2\data\RECORDING_EPHYS.xlsx');
Tall = Tall(Tall.include>0,:);

% Re-spike sort and overwrite kilosort output?
overwrite_ks = 0;

% Re-calculate GMM based on manually curated sessions, and apply to all sessions?
overwrite_gmm = 0;

% Data folder
data_dir_local = 'D:\data';  % faster to process on locally saved files
data_dir_remote = 'Z:\Hannah\Ephys\Project2'; % Location to backup to

%% Set up

% Only process data that is actually present!
fnames = dir(data_dir_local);
fnames = {fnames.name};
Tlocal = Tall(ismember(Tall.filename, fnames),:);

% Save directory function -- ks results for each session saved here
ks_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');

% Repository folder
scripts_dir = fileparts(which(mfilename));
spikesort_hp_dir = fileparts(scripts_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src')))

% Add kilosort2.0 directory (from: https://github.com/MouseLand/Kilosort/releases/tag/v2.0)
code_dir = fileparts(spikesort_hp_dir); %
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))

% Common results folder (GMM etc)
results_dir = fullfile(fileparts(scripts_dir), 'results');
addpath(results_dir);

% warning("off","parallel:gpu:device:DeviceDeprecated");

%% Backup from local SSD to server
exceptions = {'kilosort2_output','\.avi$','raw*'}; % Back this up specifically later
runmode = 2; % Don't delete anything!
for ii = 1:height(Tlocal)
    fprintf('\n%s\n',Tlocal.filename{ii})
    root_dir_local = fullfile(data_dir_local, Tlocal.filename{ii});  % e.g. D:\data\HC11_230129
    root_dir_remote = fullfile(data_dir_remote, Tlocal.filename{ii}); % e.g. Z:\Hannah\ephys\HC11_230129
    if ~exist(root_dir_local,'dir') || strcmp(root_dir_local, root_dir_remote); continue; end
    dirbackup(root_dir_local, root_dir_remote, runmode, exceptions)
end

%% Backup from server to local!
exceptions = {'kilosort2_output','\.avi$','raw*'}; % Back this up specifically later
runmode = 2; % Don't delete anything!
for ii = 1:height(Tlocal)
    fprintf('\n%s\n',Tlocal.filename{ii})
    root_dir_local = fullfile(data_dir_local, Tlocal.filename{ii});  % e.g. D:\data\HC11_230129
    root_dir_remote = fullfile(data_dir_remote, Tlocal.filename{ii}); % e.g. Z:\Hannah\ephys\HC11_230129
    if ~exist(root_dir_local,'dir') || strcmp(root_dir_local, root_dir_remote); continue; end
    dirbackup(root_dir_remote, root_dir_local, runmode, exceptions)
end

%% Run Kilosort2
for ii = 1:height(Tlocal)
    % Root dir on server (Z:\Hannah\ephys\HC11_230129 etc)
    root_dir_local = fullfile(data_dir_local, Tlocal.filename{ii});
    
    % Find the binary directory
    raw_dir_temp = dir(fullfile(root_dir_local, 'raw*'));
    if isempty(raw_dir_temp); continue; end
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    
    % Create save directory
    ks_dir = ks_dir_fun(root_dir_local);
    
    % If it exists already, either overwrite or skip
    if ~isempty(dir(fullfile(ks_dir,'*.npy')))
        if overwrite_ks
            fprintf('Overwriting %s\n', Tlocal.filename{ii})
            try
                rmdir(ks_dir,'s')
            catch
            end
        else
            fprintf('Already sorted %s\n', Tlocal.filename{ii})
            continue;
        end
    end
    mkdir(ks_dir)
    disp(Tlocal(ii,:))
    
    % set up options
    ops = [];
    ops.badChannels = eval(Tlocal.bad_chan{ii});
    ops.chanMap = [Tlocal.probe{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);
    
    % Run and save results. Note: ops saved in rez.mat
    fprintf('\n Running Kilosort on directory %s \n', raw_dir),
    runSingleKilosort(raw_dir, ks_dir, ops);
    
    % Copy kilosort output to the server
    root_dir_remote = fullfile(data_dir_remote, Tlocal.filename{ii});
    ks_dir_remote = ks_dir_fun(root_dir_remote);
    copyfile(ks_dir, ks_dir_remote)
    
    % Edit params.py to reflect new amplifier.dat location
    fid = fopen(fullfile(ks_dir,'params.py'), 'r'); % Read params.py from local copy
    a = {}; while ~feof(fid); a{end+1} = fgets(fid);  pause(.0001);  end
    fclose(fid);
    a{1} = sprintf("dat_path = '%s'\n", strrep(fullfile(root_dir_remote, raw_dir_temp.name, 'amplifier.dat'),'\','/'));
    fid = fopen(fullfile(ks_dir_remote,'params.py'), 'w'); % Write correct params.py to remote server
    for jj = 1:length(a); fprintf(fid,a{jj});      end
    fclose(fid);
    fprintf('Results copied to server!\n')
    
end

%% Go manually label results in Phy!
% In Anaconda prompt:
%   conda activate phy_env
% Navigate to folder with kilosort results (Type Z: to change to engram drive. then cd path. )
%   phy template-gui params.py

%% Analyze waveforms on local server
% analyzeBatchKilosort(Tlocal, data_dir_remote, ks_dir_fun, results_dir, overwrite_gmm)
analyzeBatchKilosort(Tlocal, data_dir_local, ks_dir_fun, results_dir, overwrite_gmm)


%% Plot spike sorting results across all sessions
plotEIClusters(Tlocal, data_dir, ks_dir_fun)

%% Backup results from local to remote
exceptions = {'params.py','phy.log','.phy'};
runmode = 2;
for ii = 1:height(Tlocal)
    fprintf('\n%s\n',Tlocal.filename{ii})
    s1 = ks_dir_fun(fullfile(data_dir_local, Tlocal.filename{ii}));
    s2 = ks_dir_fun(fullfile(data_dir_remote, Tlocal.filename{ii}));
    dirbackup(s1, s2, runmode, exceptions)
end
% Backup results from remote to local 
exceptions = {'params.py','phy.log','.phy'};
runmode = 2;
for ii = 1:height(Tlocal)
    fprintf('\n%s\n',Tlocal.filename{ii})
    s1 = ks_dir_fun(fullfile(data_dir_local, Tlocal.filename{ii}));
    s2 = ks_dir_fun(fullfile(data_dir_remote, Tlocal.filename{ii}));
    dirbackup(s2, s1, runmode, exceptions)
end


%% Optional plotting of # of cells and E/I ratio
plotBasicResults(Tlocal, data_dir_remote, ks_dir_fun)
