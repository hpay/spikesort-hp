function info = readOpenEphysInfo(rec_dir, stream_name)
% READOPENEPHYSINFO  Read the structure.oebin of an openEphys binary-format recording
%
% info = readOpenEphysInfo(rec_dir, stream_name)
%   rec_dir: ...\Record Node XXX\experimentN\recordingN
%   stream_name: continuous folder name (default 'Acquisition_Board-100.acquisition_board')
%
% info fields: fbinary, nCh, fs, nSamp, channel_names, bit_volts, units,
%   ephys_idx (1-indexed CH* channels), adc_idx (1-indexed ADC* channels)

if ~exist('stream_name','var') || isempty(stream_name)
    stream_name = 'Acquisition_Board-100.acquisition_board';
end

oebin = jsondecode(fileread(fullfile(rec_dir, 'structure.oebin')));
folders = cellfun(@(f) strip(f, 'right', '/'), {oebin.continuous.folder_name}, 'Uni', 0);
i_stream = find(strcmp(folders, stream_name));
if isempty(i_stream); error('Stream %s not found in structure.oebin', stream_name); end
stream = oebin.continuous(i_stream);

info.fbinary = fullfile(rec_dir, 'continuous', stream_name, 'continuous.dat');
info.nCh = stream.num_channels;
info.fs = stream.sample_rate;
info.channel_names = {stream.channels.channel_name}';
info.bit_volts = [stream.channels.bit_volts]';
info.units = {stream.channels.units}';
info.ephys_idx = find(startsWith(info.channel_names, 'CH'))';
info.adc_idx = find(startsWith(info.channel_names, 'ADC'))';

f = dir(info.fbinary);
info.nSamp = f.bytes / (2*info.nCh);
if info.nSamp ~= round(info.nSamp); error('File size is not a whole number of samples'); end
