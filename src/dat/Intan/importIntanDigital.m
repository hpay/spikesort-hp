function datout = importIntanDigital(filepath, option_raw_dig)
% Converts Intan digital input to continuous data structure
%
% Input:
%    filename
% 
% Output:
%   datout - dat object with fields:
%     'data' - column vector of numeric data
%     'chanlabels' - string labels for each channel
%     'chanvals' - numeric values associated with each
%         channel
%     'samplerate' - derived, not claimed (samples/second, double-precision)
%     'tstart','tend' - time of first/last sample; in seconds.
%     'units' - descriptive string giving units of data.

%
% Dependencies: readIntanPerType
%
% Hannah Payne 2019

if ~exist('filepath','var')
    [filepath] = uigetdir('Select a folder with single file per data type data');
    if (filepath == 0);  return; end

end

if ~exist('option_raw_dig','var')
    option_raw_dig = 0;
end

% Read raw data and header
D = [];
option_digital_only = 1;
[D.amplifier_data, D.header, D.board_dig_in_data] = readIntanPerType(filepath, [], [], [], option_digital_only);
samplerate = D.header.sample_rate;
nt = size(D.board_dig_in_data,2);

% Assign time stamps
D.t_dig = (0:nt-1)/samplerate;

% Read timestamps
% fid = fopen(fullfile(filepath,'time.dat'));
% D.t_dig = fread(fid,'int32');
% fclose(fid);
% if any(diff(D.t_dig)~=1)
%     warning('Check timestampes')
% end
% D.t_dig = D.t_dig;

%% Save the board digital input data in dat structure
if isfield(D.header,'board_dig_in_channels') && ~isempty(D.header.board_dig_in_channels)
    iunits = 's'; %
    
    % Get time and unit info
    itstart = 0; % (s) Time since start of recording
    itend = itstart + (size(D.board_dig_in_data,2)-1)/samplerate; 
    
    nchans_dig = length(D.header.board_dig_in_channels);
    
    % Preallocate
    datout = dat(nchans_dig);

    for jj = 1:nchans_dig
        
        ichanlabel = D.header.board_dig_in_channels(jj).custom_channel_name;
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
            
            datout(jj) = dat(iData,ichanlabel,ichanval,'dig',itstart,itend,iunits);
        else
            iData = D.board_dig_in_data(jj,:)>0;
            datout(jj) = dat(iData,ichanlabel,ichanval,samplerate,itstart,itend,iunits);
            
        end
    end
end


%% Create info field with all other information
datout(1).info = D.header;

