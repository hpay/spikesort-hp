function rawdata = readIntanPerChannel(filename, duration_ind, start_inds, ndownsample)
% readIntanPerChannel Returns raw voltage trace in uV
%
% rawdata = readIntanPerChannel(filename)
% rawdata = readIntanPerChannel(filename, duration_ind, start_inds)
%
% Example: rawdata = readIntanPerChannel('Z:\Hannah\ephys\HT09_181221\raw_181221_112114\amp-A-033.dat')

fid = fopen(filename,'r');
if ~exist('ndownsample','var')
nskip = 0;
else
nskip = ndownsample-1;
end
if nargin==1 || isempty(duration_ind)
    rawdata = single(fread(fid, 'int16', nskip*2))*0.195;

elseif nargin==2 || isempty(start_inds) 
    rawdata = single(fread(fid, duration_ind, 'int16', nskip*2))*0.195;

elseif nargin>=3
    rawdata = NaN(duration_ind, length(start_inds));
    for ii = 1:length(start_inds)
        status = fseek(fid,  2*start_inds(ii), -1);
        if status
            error('fseek not successful in readIntanPerChannel, check starting time')
        end
        currdata = fread(fid, duration_ind, 'int16', nskip*2)*0.195;
        if length(currdata) == duration_ind  
        rawdata(:,ii) = currdata;
        end
    end
end
fclose(fid);