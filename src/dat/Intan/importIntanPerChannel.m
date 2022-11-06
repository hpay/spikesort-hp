function [D, info] = importIntanPerChannel(foldername, tduration, tstart, ch_num, printoutput)
%  [D, h] = importIntanPerChannel(foldername, tduration, tstart, ch_num)
% Load the raw data from all single data channels (one file per channel
% format) present in a given folder
% Returns data in D in dat format, info in h


if exist('ch_num','var') && ~isempty(ch_num)
    amp_files = dir(fullfile(foldername,sprintf('amp*%03i.dat', ch_num)));
else
    amp_files = dir(fullfile(foldername,'amp*.dat'));
end
if ~exist('printoutput','var')
    printoutput = 1;
end

try
    info = readIntanInfo(fullfile(foldername,'info.rhd'),printoutput);
    samplerate = info.sample_rate;
catch
    try
        info = load(fullfile(foldername,'info.mat'));
        samplerate = info.header.sample_rate;
    catch
        samplerate = 1;
    end
    
end
units = 'uV';

duration_ind = [];
start_ind = [];
if exist('tduration','var') && ~isempty(tduration)
    duration_ind = round(tduration*samplerate);
end
if exist('tstart','var') && ~isempty(tstart)
    start_ind = round(tstart*samplerate);
else
    tstart = 0;
end

if printoutput
    fprintf('Loading channels ')
end

for ii = 1:length(amp_files)
    if printoutput
        fprintf('%i %s. ', (ii),amp_files(ii).name)
    end
    rawdata = readIntanPerChannel(fullfile(foldername, amp_files(ii).name), duration_ind, start_ind);
    chanlabel = amp_files(ii).name(5:9);
    D(ii) = dat(rawdata, chanlabel, ii, samplerate, tstart, [], units);
    
end

if printoutput
    fprintf('\n')
end