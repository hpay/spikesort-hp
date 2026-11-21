function datout = importIntanDigital(filepath, option_raw_dig)
% Converts Intan digital input to continuous data structure
%
% Input:
%    filename
% 
% Output:
%   datout - struct array, one element per channel, with fields:
%     'data' - N x 2 [rise fall] times (s), one pulse per row; first rise is NaN if the
%              line starts high, last fall is NaN if it ends high (raw 0/1 trace if option_raw_dig)
%     'chanlabel' - channel name
%     'chanval' - channel number
%     'samplerate' - 'dig' for edge times, else samples/s
%     'tstart','tend' - time of first/last sample (s)
%     'units' - 's'
%     'info' - (first element only) file header

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
    
    datout = struct([]);
    for jj = 1:nchans_dig
        ichanlabel = D.header.board_dig_in_channels(jj).custom_channel_name;
        ichanval = D.header.board_dig_in_channels(jj).native_order;
        if ~option_raw_dig
            x = D.board_dig_in_data(jj,:);
            rise = D.t_dig(diff(x)>0)'; fall = D.t_dig(diff(x)<0)';
            if ~isempty(fall) && (isempty(rise) || fall(1) < rise(1)); rise = [NaN; rise]; end % starts high
            if numel(fall) < numel(rise); fall(end+1,1) = NaN; end % ends high
            iData = [rise fall];
            isamplerate = 'dig';
        else
            iData = (D.board_dig_in_data(jj,:)>0)';
            isamplerate = samplerate;
        end
        datout(jj).data = iData;
        datout(jj).chanlabel = ichanlabel;
        datout(jj).chanval = ichanval;
        datout(jj).samplerate = isamplerate;
        datout(jj).tstart = itstart;
        datout(jj).tend = itend;
        datout(jj).units = iunits;
    end
end


%% Create info field with all other information
datout(1).info = D.header;

