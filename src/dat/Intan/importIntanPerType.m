function datout = importIntanPerType(filepath, tduration, tstart,  notch_filter_frequency, read_dig, option_raw_dig, option_type)
% Converts Intan single file per channel type to continuous data structure
%
% Input:
%    filename
% 	 (tduration): (s)  how many seconds of data to read
%    (tstart): (s) when to start readibng
%    (notch_filter_frequency): (Hz) notch filter e.g. 60 Hz
%    (read_dig): include digital
%    (option_raw_dig): don't convert digital to events, keep it as time
% Output:
%   datout - dat object with fields:
%     'data' - column vector of numeric data
%     'chanlabels' - string labels for each channel
%     'chanvals' - numeric values associated with each
%         channel
%     'samplerate' - derived, not claimed (samples/second, double-precision)
%     'tstart','tend' - time of first/last sample; in seconds.
%     'units' - descriptive string giving units of data.
%     ('wavemark') - spike identity
%     ('waveform') for waveform spk channel from spike2
%
%
% Dependencies: readIntanSingleFile
%
% Hannah Payne 2019
%
% TODO: add notch filter option
% TODO: add non-amplifier channels!
% TODO: add option_Raw_dig

if ~exist('filepath','var')
    [filepath] = uigetdir('Select a folder with single file per data type data');
    if (filepath == 0);  return; end

end

if ~exist('tduration','var')
    tduration = []; 
end

if ~exist('tstart','var')
    tstart = []; 
end

if ~exist('notch_filter_frequency','var')
    notch_filter_frequency = 0;
end

read_dig = 1; % read digital data as well or not

if ~exist('option_raw_dig','var')
    option_raw_dig = 0;
end

% Read raw data and header
D = [];
[D.amplifier_data, D.header] = readIntanPerType(filepath, tduration, tstart, notch_filter_frequency, read_dig);


% Read timestamps
fid = fopen(fullfile(filepath,'time.dat'));
if ~isempty(tstart)
    ind_start = round(tstart*D.header.sample_rate);
    fseek(fid, ind_start*4, -1);
end
if ~isempty(tduration)
    ind_duration= round(tduration*D.header.sample_rate);
    D.t_amplifier = fread(fid,ind_duration,'int32');
else
    D.t_amplifier = fread(fid,'int32');
end
fclose(fid);
if any(diff(D.t_amplifier)~=1)
    warning('Check timestampes')
end

% Get basic facts
nchans = length(D.header.amplifier_channels);
nsamples = length(D.t_amplifier);

curr_chan = 1;


% Preallocate object array
if nchans>0
    datout(nchans) = dat;
end

% To prevent out of memory errors
if isfield(D,'amplifier_data')
    D.amplifier_data = single(D.amplifier_data);
end
if isfield(D,'board_dig_in_data')
    D.board_dig_in_data = D.board_dig_in_data>0;
end

%% Store amplifier data in dat structure
iunits = 'uV'; % Microvolts
samplerate = D.header.sample_rate;

% Get time and unit info
% itstart = D.t_amplifier(1); % (s) Time since start of recording
itstart = 0;
itend = itstart + (nsamples-1)/samplerate;
for ii = 1:nchans
    
    ichanlabel = D.header.amplifier_channels(ii).native_channel_name;
    ichanval = D.header.amplifier_channels(ii).native_order;
    
    % Convert current data to single
    iData = D.amplifier_data(ii,:);
    
    % Make the dat structure channel
    datout(curr_chan) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
    curr_chan = curr_chan+1;
end



%% Save the board digital input data in dat structure
%%{
nchans_dig = 0;
if isfield(D.header,'board_dig_in_channels') && ~isempty(D.header.board_dig_in_channels)
    iunits = 's'; %
    
    % Get time and unit info
    itstart = D.t_dig(1); % (s) Time since start of recording
    itend = itstart + (size(D.board_dig_in_data,2)-1)/samplerate; % == D.t_dig(end)
    
    nchans_dig = length(D.header.board_dig_in_channels);
    
    for jj = 1:nchans_dig
        
        ichanlabel = D.header.board_dig_in_channels(jj).native_channel_name;
        ichanval = D.header.board_dig_in_channels(jj).native_order;
        
        if ~option_raw_dig
            
            iData_rise = D.t_dig(diff(D.board_dig_in_data(jj,:))>0);
            iData_fall = D.t_dig(diff(D.board_dig_in_data(jj,:))<0);
            
            iData_rise = iData_rise(:);
            iData_fall = iData_fall(:);
            
            % First column is rise times, second column is fall times
            try
                if iData_rise(1)>iData_fall(1)
                    iData_rise = [NaN; iData_rise];
                end
                if iData_fall(end)<iData_rise(end)
                    iData_fall = [iData_fall; NaN];
                end
            end
            iData = [iData_rise iData_fall];
            
            datout(curr_chan) = dat(iData,ichanlabel,ichanval,'dig',itstart,itend,iunits);
        else
            iData = D.board_dig_in_data(jj,:)>0;
            datout(curr_chan) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
            
        end
         curr_chan = curr_chan+1;

    end
end

%% Save the board adc input data in dat structure
nchans_board_adc = 0;
if isfield(D.header,'board_adc_channels') && ~isempty(D.header.board_adc_channels)
    iunits = 'V';
    samplerate = D.frequency_parameters.board_adc_sample_rate;
    
    nchans_board_adc = length(D.header.board_adc_channels);
    
    for jj = 1:nchans_board_adc
        
        ichanlabel = D.header.board_adc_channels(jj).native_channel_name;
        ichanval = D.header.board_adc_channels(jj).native_order;
        
        iData = D.board_adc_data(jj,:);
        
        datout(curr_chan) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
        curr_chan = curr_chan+1;

    end
end

%% Save the aux input data in dat format
if isfield(D.header,'aux_input_channels') && ~isempty(D.header.aux_input_channels)
    iunits = 'V';
    samplerate = D.frequency_parameters.aux_input_sample_rate;
    
    nchans_board_aux = length(D.header.aux_input_channels);
    
    for jj = 1:nchans_board_aux
        
        ichanlabel = D.header.aux_input_channels(jj).native_channel_name;
        ichanval = D.header.aux_input_channels(jj).native_order;
        
        iData = D.aux_input_data(jj,:);
        
        datout(curr_chan) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
        curr_chan = curr_chan+1;

    end
end

%}

%% Create info field with all other information
datout(1).info = D.header;

