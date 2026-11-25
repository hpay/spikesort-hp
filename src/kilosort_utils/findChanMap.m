function chan_map_file = findChanMap(ks_dir)
% FINDCHANMAP  Channel map file (in src/kilosort_config_files) a finished Kilosort 2 or 4 sort used: the one
% whose x, y positions match the sort's channel_positions.npy for the channels it kept (channel_map.npy,
% 0-indexed raw channel numbers). Returns '' if no map matches.

ch = double(readNPY(fullfile(ks_dir, 'channel_map.npy'))) + 1;
pos = double(readNPY(fullfile(ks_dir, 'channel_positions.npy')));
map_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'kilosort_config_files');
files = dir(fullfile(map_dir, '*.mat'));
match = {};
for f = files'
    m = load(fullfile(f.folder, f.name));
    if ~all(isfield(m, {'chanMap', 'xcoords', 'ycoords'})); continue; end
    [found, loc] = ismember(ch(:), m.chanMap(:));
    if all(found) && max(abs([m.xcoords(loc(:)), m.ycoords(loc(:))] - pos), [], 'all') < 1e-3
        match{end + 1} = f.name; %#ok<AGROW>
    end
end
chan_map_file = '';
if ~isempty(match); chan_map_file = match{1}; end
if numel(match) > 1; warning('%s matches several channel maps: %s; using %s', ks_dir, strjoin(match, ', '), match{1}); end
end
