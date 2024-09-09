function plotEIClusters(T, data_dir, ks_dir_fun)
%
% Inputs:
%   T - A table containing session metadata (e.g., filenames, inclusion criteria).
%   data_dir - A string specifying the root directory where data files are stored.
%   ks_dir_fun - A function handle that constructs the path to the Kilosort output directory for each session.
%
% Outputs:
%   Generates a plot displaying cluster analysis results, including waveform data and GMM clustering results.
%
% Example usage
% Tall = readtable('C:\Users\Hannah\Dropbox\alab\Code\project2\data\RECORDING_EPHYS.xlsx');
% Tall = Tall(Tall.include>0,:);
% data_dir_local = 'D:\data';  % faster to process on locally saved files
% ks_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');


stats_all = [];
wv_all = [];
id_all = [];
good_all = [];
session_all = [];
channel_all = [];
unitID_all = [];
for ii = height(T):-1:1
    ks_dir = ks_dir_fun(fullfile(data_dir, T.filename{ii}));
    gmm_curr = load(fullfile(ks_dir,'gmm_result.mat'),'id','labels','stats');
    wvstruct_curr = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')),'wvStruct');
    stats_all = [stats_all; gmm_curr.stats];
    wv_all = [wv_all; wvstruct_curr.mxWF];
    id_all = [id_all; wvstruct_curr.typeLabels];
    good_all = [good_all; wvstruct_curr.goodLabels];
    session_all = [session_all; repmat(T.filename(ii), length(gmm_curr.id),1)];
    channel_all = [channel_all; wvstruct_curr.max_site];
    unitID_all= [unitID_all; wvstruct_curr.goodIDs];
end
nt = size(wvstruct_curr.mxWF,2);
dt = wvstruct_curr.spkDur/nt;
tWF = -wvstruct_curr.spkOffset + (0:nt-1)*dt;
Ranatomy = getAnatomy(session_all, channel_all, T); % Store AP/ML/DV info
 
plotGMMresults(tWF, wv_all, stats_all, id_all, good_all, session_all, Ranatomy, channel_all, unitID_all)
drawnow
