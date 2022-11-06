function [datout, tamplifier] = importIntanRHD(filename, notch_filter_frequency, option_raw_dig)
%% Converts .rhd to continuous data structure
%
% Input:
%    filename
%    notch_filter_frequency (optional): 60 Hz
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
% digout - dat object with digital data
% adcout - dat object with analog data
%
% Dependencies: readIntanRHD
%
% Hannah Payne April 2018

if ~exist('notch_filter_frequency','var')
    notch_filter_frequency = 0;
end

if ~exist('option_raw_dig','var')
    option_raw_dig = 0;
end

D = readIntanRHD(filename, notch_filter_frequency);
tamplifier = D.t_amplifier;

% Get basic facts
nchans = length(D.header.amplifier_channels);

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
samplerate = D.frequency_parameters.amplifier_sample_rate;

% Get time and unit info
itstart = D.t_amplifier(1); % (s) Time since start of recording
itend = itstart + (size(D.amplifier_data,2)-1)/samplerate;

for ii = 1:nchans
    
    ichanlabel = D.header.amplifier_channels(ii).native_channel_name;
    ichanval = D.header.amplifier_channels(ii).native_order;
    
    % Convert current data to single
    iData = D.amplifier_data(ii,:);
    
    % Make the dat structure channel
    datout(ii) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
    
end

%% Save the board digital input data in dat structure
nchans_dig = 0;
if isfield(D.header,'board_dig_in_channels') && ~isempty(D.header.board_dig_in_channels)
    iunits = 's'; %
    samplerate = D.frequency_parameters.board_dig_in_sample_rate;
    
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
            
            datout(nchans+jj) = dat(iData,ichanlabel,ichanval,'dig',itstart,itend,iunits);
        else
            iData = D.board_dig_in_data(jj,:)>0;
            datout(nchans+jj) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
            
        end
        
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
        
        datout(nchans + nchans_dig + jj) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
        
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
        
        datout(nchans + nchans_dig + nchans_board_adc + jj) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
        
    end
end

%% Create info field with all other information
datout(1).info = D.header;

