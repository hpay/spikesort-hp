function DIN = importOpenEphysDigital(fbinary, nCh, chan_idx, thresh, fs, channel_names)
% IMPORTOPENEPHYSDIGITAL  Threshold channels of an interleaved int16 binary (e.g. openEphys
% ADC inputs carrying TTLs) and return rise/fall times in the same struct format as
% importIntanDigital (DIN.mat). Reads in chunks; only reads fbinary.
%
% DIN = importOpenEphysDigital(fbinary, nCh, chan_idx, thresh, fs, channel_names)
%   chan_idx: 1-indexed channels in the file
%   thresh: threshold in raw int16 counts (high = x > thresh)
%   DIN(k).data: N x 2 [rise fall] times (s), NaN-padded like importIntanDigital. As in
%       importIntanDigital, an edge's time is that of the last sample before it,
%       t = (0-indexed sample number)/fs
%   DIN(1).info: fbinary, thresh, initial_state (state of each channel at the first sample)

chunk_samples = fs*60; % 1 min per read
n_chan = numel(chan_idx);
f = dir(fbinary);
n_samp = f.bytes / (2*nCh);

fid = fopen(fbinary, 'r');
cleanup = onCleanup(@() fclose(fid));
rise = cell(n_chan, 1);
fall = cell(n_chan, 1);
n_done = 0;
prev = [];
t = tic;
while n_done < n_samp
    n_read = min(chunk_samples, n_samp - n_done);
    x = fread(fid, [nCh n_read], '*int16');
    b = x(chan_idx,:) > thresh;
    if isempty(prev); initial_state = b(:,1); prev = b(:,1); end
    d = diff(int8([prev b]), 1, 2); % column j = transition out of 0-indexed sample n_done + j - 2
    for k = 1:n_chan
        rise{k} = [rise{k}; (find(d(k,:) > 0)' + n_done - 2) / fs];
        fall{k} = [fall{k}; (find(d(k,:) < 0)' + n_done - 2) / fs];
    end
    prev = b(:,end);
    n_done = n_done + n_read;
    if mod(n_done, 10*chunk_samples) == 0 || n_done == n_samp
        fprintf('Thresholded %.1f / %.1f min (%.0f s)\n', n_done/fs/60, n_samp/fs/60, toc(t));
    end
end
DIN = struct([]);
for k = 1:n_chan
    r = rise{k}; f = fall{k};
    if ~isempty(f) && (isempty(r) || f(1) < r(1)); r = [NaN; r]; end % starts high
    if numel(f) < numel(r); f(end+1,1) = NaN; end % ends high
    DIN(k).data = [r f];
    DIN(k).chanlabel = channel_names{k};
    DIN(k).chanval = chan_idx(k);
    DIN(k).samplerate = 'dig';
    DIN(k).tstart = 0;
    DIN(k).tend = (n_samp - 1)/fs;
    DIN(k).units = 's';
end
DIN(1).info = struct('fbinary', fbinary, 'thresh', thresh, 'initial_state', initial_state);
