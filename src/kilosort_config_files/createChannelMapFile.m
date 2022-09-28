%  create a channel map file

% H6 probe
T = readtable('D:\hannah\Dropbox\alab\Analysis\SILICON PROBE MAP H6.xlsx');

T = sortrows(sortrows(T,'ShankPosition'),'ShankLetter'); % Sort by Intan #
% T([8 12 14 16 33:end],:) = []; %***Optionally remove unrecorded channels
% T([33:end],:) = []; %***Optionally remove unrecorded channels
Nchannels = height(T);


T = sortrows(T,'IntanChannel'); % Sort by Intan #


xcoords   = T.xpos;
ycoords   = T.ypos;
kcoords   = cellfun(@(x) x-'A', T.ShankLetter)+1; % grouping of channels (i.e. tetrode groups)

fs = 30000; % sampling frequency

connected = true(Nchannels, 1); % Change if needed
xcoords = xcoords(connected);
ycoords = ycoords(connected);
kcoords = kcoords(connected);
connected = connected(connected);

chanMap   = 1:Nchannels;
% chanMap   = T.IntanChannel+1;
chanMap0ind = chanMap - 1;

% chanMap = chanMap(connected);
% chanMap0ind = chanMap0ind(connected);

save('D:\hannah\Dropbox\alab\Code\Kilosort-2.0\configFiles\chanMap_H6.mat', ... % _HC05_220811
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%% Original examples
Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 25000; % sampling frequency
save('C:\DATA\Spikes\20150601_chan32_4_900s\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%%

Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 25000; % sampling frequency

save('C:\DATA\Spikes\Piroska\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 