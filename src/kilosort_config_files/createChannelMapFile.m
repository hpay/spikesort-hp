function createChannelMapFile(xlsx_file, out_file)
% CREATECHANNELMAPFILE  Make a Kilosort channel map .mat from a probe map spreadsheet
% (not needed for spikeinterface, just for running KS2 directly in matlab)
%
% Spreadsheet needs columns IntanChannel (0-indexed channel in the raw .dat),
% ShankLetter, ShankPosition, xpos, ypos. If it also has shankPosition_1
% (0-indexed virtual shank, e.g. front/back of the double-sided F probe treated
% as separate shanks offset in x), that is used for kcoords instead of ShankLetter.
%
% Examples:
%   createChannelMapFile('SILICON PROBE MAP H6.xlsx', 'H6.mat')
%   createChannelMapFile('Z:\Meghana\probe maps\SILICON PROBE MAP F-D_12Shank1200Microns.xlsx', 'F_D_12shank.mat')

T = readtable(xlsx_file, 'VariableNamingRule', 'modify');
has_virtual_shank = ismember('shankPosition_1', T.Properties.VariableNames);

% Sort by shank & position
if has_virtual_shank
    T = sortrows(T, {'shankPosition_1', 'ShankPosition'});
    kcoords = T.shankPosition_1 + 1;
else
    T = sortrows(sortrows(T, 'ShankPosition'), 'ShankLetter');
    kcoords = cellfun(@(x) x-'A', T.ShankLetter) + 1; % grouping of channels (i.e. tetrode groups)
end
n_channels = height(T);

xcoords = T.xpos;
ycoords = T.ypos;
fs = 30000; % sampling frequency
connected = true(n_channels, 1); % Change if needed - but better to mask bad channels later
chanMap = T.IntanChannel + 1;
chanMap0ind = chanMap - 1;

save(out_file, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
if ismember('Back_Front', T.Properties.VariableNames) % double-sided probe: physical shank letter and side
    shank_letter = char(T.ShankLetter); side = char(T.Back_Front);
    save(out_file, 'shank_letter', 'side', '-append')
end
fprintf('Saved %d-channel map to %s\n', n_channels, out_file)
