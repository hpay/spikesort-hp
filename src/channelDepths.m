function [depth, shank, tip_depth] = channelDepths(session, raw_chan, data_dir, chan_map_file)
% CHANNELDEPTHS  Depth below the brain surface and shank of raw channels (CH#, 1-indexed) in a
% session: probe tip position in the brain (DataNotebook-MH.xlsx, Ephys sheet, "Probe position
% in the brain") minus the channel's height above the tip (ycoords in the channel map).
% Channels not in the map (e.g. NaN) get NaN.

m = load(chan_map_file);
T_ephys = readtable(fullfile(data_dir, 'DataNotebook-MH.xlsx'), 'Sheet', 'Ephys', 'VariableNamingRule', 'preserve');
tip_depth = T_ephys{strcmp(T_ephys.Folder, session), 'Probe position in the brain = absolute - 250 uM'};
[in_map, map_idx] = ismember(raw_chan(:), m.chanMap);
depth = nan(size(map_idx)); shank = nan(size(map_idx));
depth(in_map) = tip_depth - m.ycoords(map_idx(in_map));
shank(in_map) = m.kcoords(map_idx(in_map));
end
