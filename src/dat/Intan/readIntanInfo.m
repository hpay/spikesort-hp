function h = readIntanInfo(filename,printoutput)
% Returns header information for an info.rhd file


if ~exist('filename','var')
    [file, path] = ...
        uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'off');
    if (file == 0);  return; end
else
    [path, file, ext] = fileparts(filename);
    if isempty(ext); ext = '.rhd'; end
    file = [file ext];
end

% Set a software 50/60 Hz notch filter with input notch_filter_frequency
if ~exist('notch_filter_frequency','var')
    notch_filter_frequency = 0;
end

% Open the file
tic;
filename = fullfile(path,file);
if ~exist('printoutput','var') || printoutput
fprintf('\nOpening file %s\n',filename)
end
fid = fopen(filename, 'r');

if fid<1; warning('Failed to open info file'); return; end

s = dir(filename);
filesize = s.bytes;

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
h.magic_number = fread(fid, 1, 'uint32');
if h.magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
h.data_file_main_version_number = fread(fid, 1, 'int16');
h.data_file_secondary_version_number = fread(fid, 1, 'int16');

if ~exist('printoutput','var') || printoutput
fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
    h.data_file_main_version_number, h.data_file_secondary_version_number);
end

if (h.data_file_main_version_number == 1)
    num_samples_per_data_block = 60;
else
    num_samples_per_data_block = 128;
end

% Read information of sampling rate and amplifier frequency settings.
h.sample_rate = fread(fid, 1, 'single');
h.dsp_enabled = fread(fid, 1, 'int16');
h.actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
h.actual_lower_bandwidth = fread(fid, 1, 'single');
h.actual_upper_bandwidth = fread(fid, 1, 'single');

h.desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
h.desired_lower_bandwidth = fread(fid, 1, 'single');
h.desired_upper_bandwidth = fread(fid, 1, 'single');

% This tells us if a software 50/60 Hz notch filter was enabled during
% the data acquisition.
h.notch_filter_mode = fread(fid, 1, 'int16');

h.desired_impedance_test_frequency = fread(fid, 1, 'single');
h.actual_impedance_test_frequency = fread(fid, 1, 'single');

% Place notes in data strucure
h.notes = struct( ...
    'note1', fread_QString(fid), ...
    'note2', fread_QString(fid), ...
    'note3', fread_QString(fid) );
    
% If data file is from GUI v1.1 or later, see if temperature sensor data
% was saved.
h.num_temp_sensor_channels = 0;
if ((h.data_file_main_version_number == 1 && h.data_file_secondary_version_number >= 1) ...
    || (h.data_file_main_version_number > 1))
    h.num_temp_sensor_channels = fread(fid, 1, 'int16');
end

% If data file is from GUI v1.3 or later, load eval board mode.
h.eval_board_mode = 0;
if ((h.data_file_main_version_number == 1 && h.data_file_secondary_version_number >= 3) ...
    || (h.data_file_main_version_number > 1))
    h.eval_board_mode = fread(fid, 1, 'int16');
end

% If data file is from v2.0 or later (Intan Recording Controller),
% load name of digital reference channel.
if (h.data_file_main_version_number > 1)
    h.reference_channel = fread_QString(fid);
end

% Place frequency-related information in data structure.
frequency_parameters = struct( ...
    'amplifier_sample_rate', h.sample_rate, ...
    'aux_input_sample_rate', h.sample_rate / 4, ...
    'supply_voltage_sample_rate', h.sample_rate / num_samples_per_data_block, ...
    'board_adc_sample_rate', h.sample_rate, ...
    'board_dig_in_sample_rate', h.sample_rate, ...
    'desired_dsp_cutoff_frequency', h.desired_dsp_cutoff_frequency, ...
    'actual_dsp_cutoff_frequency', h.actual_dsp_cutoff_frequency, ...
    'dsp_enabled', h.dsp_enabled, ...
    'desired_lower_bandwidth', h.desired_lower_bandwidth, ...
    'actual_lower_bandwidth', h.actual_lower_bandwidth, ...
    'desired_upper_bandwidth', h.desired_upper_bandwidth, ...
    'actual_upper_bandwidth', h.actual_upper_bandwidth, ...
    'notch_filter_frequency', notch_filter_frequency, ...
    'desired_impedance_test_frequency', h.desired_impedance_test_frequency, ...
    'actual_impedance_test_frequency', h.actual_impedance_test_frequency );

% Define data structure for spike trigger settings.
spike_trigger_struct = struct( ...
    'voltage_trigger_mode', {}, ...
    'voltage_threshold', {}, ...
    'digital_trigger_channel', {}, ...
    'digital_edge_polarity', {} );

new_trigger_channel = struct(spike_trigger_struct);
h.spike_triggers = struct(spike_trigger_struct);

% Define data structure for data channels.
channel_struct = struct( ...
    'native_channel_name', {}, ...
    'custom_channel_name', {}, ...
    'native_order', {}, ...
    'custom_order', {}, ...
    'board_stream', {}, ...
    'chip_channel', {}, ...
    'port_name', {}, ...
    'port_prefix', {}, ...
    'port_number', {}, ...
    'electrode_impedance_magnitude', {}, ...
    'electrode_impedance_phase', {} );

new_channel = struct(channel_struct);

% Create structure arrays for each type of data channel.
h.amplifier_channels = struct(channel_struct);
h.aux_input_channels = struct(channel_struct);
h.supply_voltage_channels = struct(channel_struct);
h.board_adc_channels = struct(channel_struct);
h.board_dig_in_channels = struct(channel_struct);
h.board_dig_out_channels = struct(channel_struct);

amplifier_index = 1;
aux_input_index = 1;
supply_voltage_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;
board_dig_out_index = 1;

% Read signal summary from data file header.

h.number_of_signal_groups = fread(fid, 1, 'int16');

for signal_group = 1:h.number_of_signal_groups
    signal_group_name = fread_QString(fid);
    signal_group_prefix = fread_QString(fid);
    signal_group_enabled = fread(fid, 1, 'int16');
    signal_group_num_channels = fread(fid, 1, 'int16');
    signal_group_num_amp_channels = fread(fid, 1, 'int16');

    if (signal_group_num_channels > 0 && signal_group_enabled > 0)
        new_channel(1).port_name = signal_group_name;
        new_channel(1).port_prefix = signal_group_prefix;
        new_channel(1).port_number = signal_group;
        for signal_channel = 1:signal_group_num_channels
            new_channel(1).native_channel_name = fread_QString(fid);
            new_channel(1).custom_channel_name = fread_QString(fid);
            new_channel(1).native_order = fread(fid, 1, 'int16');
            new_channel(1).custom_order = fread(fid, 1, 'int16');
            signal_type = fread(fid, 1, 'int16');
            channel_enabled = fread(fid, 1, 'int16');
            new_channel(1).chip_channel = fread(fid, 1, 'int16');
            new_channel(1).board_stream = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
            new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
            new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
            new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
            new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
            
            if (channel_enabled)
                switch (signal_type)
                    case 0
                        h.amplifier_channels(amplifier_index) = new_channel;
                        h.spike_triggers(amplifier_index) = new_trigger_channel;
                        amplifier_index = amplifier_index + 1;
                    case 1
                        h.aux_input_channels(aux_input_index) = new_channel;
                        aux_input_index = aux_input_index + 1;
                    case 2
                        h.supply_voltage_channels(supply_voltage_index) = new_channel;
                        supply_voltage_index = supply_voltage_index + 1;
                    case 3
                        h.board_adc_channels(board_adc_index) = new_channel;
                        board_adc_index = board_adc_index + 1;
                    case 4
                        h.board_dig_in_channels(board_dig_in_index) = new_channel;
                        board_dig_in_index = board_dig_in_index + 1;
                    case 5
                        h.board_dig_out_channels(board_dig_out_index) = new_channel;
                        board_dig_out_index = board_dig_out_index + 1;
                    otherwise
                        error('Unknown channel type');
                end
            end
            
        end
    end
end

% Summarize contents of data file.
num_h.amplifier_channels = amplifier_index - 1;
num_h.aux_input_channels = aux_input_index - 1;
num_h.supply_voltage_channels = supply_voltage_index - 1;
num_h.board_adc_channels = board_adc_index - 1;
num_h.board_dig_in_channels = board_dig_in_index - 1;
num_h.board_dig_out_channels = board_dig_out_index - 1;
h.num_h = num_h;

if ~exist('printoutput','var') || printoutput
fprintf(1, 'Found %d amplifier channel%s.\n', ...
    num_h.amplifier_channels, plural(num_h.amplifier_channels));
fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
    num_h.aux_input_channels, plural(num_h.aux_input_channels));
fprintf(1, 'Found %d supply voltage channel%s.\n', ...
    num_h.supply_voltage_channels, plural(num_h.supply_voltage_channels));
fprintf(1, 'Found %d board ADC channel%s.\n', ...
    num_h.board_adc_channels, plural(num_h.board_adc_channels));
fprintf(1, 'Found %d board digital input channel%s.\n', ...
    num_h.board_dig_in_channels, plural(num_h.board_dig_in_channels));
fprintf(1, 'Found %d board digital output channel%s.\n', ...
    num_h.board_dig_out_channels, plural(num_h.board_dig_out_channels));
fprintf(1, 'Found %d temperature sensor channel%s.\n', ...
    h.num_temp_sensor_channels, plural(h.num_temp_sensor_channels));
fprintf(1, '\n');
end

fclose(fid);

return

function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return


function s = plural(n)

% s = plural(n)
% 
% Utility function to optionally plurailze words based on the value
% of n.

if (n == 1)
    s = '';
else
    s = 's';
end

return
