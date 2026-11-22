function [S, f] = channelSpectra(fbinary, nCh, fs, chan_idx, bit_volts, opts)
% CHANNELSPECTRA  Multitaper power spectrum (Chronux mtspectrumc) of each channel of an
% interleaved int16 binary (read only). Spectra are computed on opts.n_seg segments of
% opts.seg_dur s spread evenly through the recording, each downsampled to opts.fs_lfp, and
% averaged over segments.
%   chan_idx: 1-indexed channels in the file
%   S: frequencies x channels power (uV^2/Hz), columns in chan_idx order; f: frequencies (Hz)
%   opts (all optional): n_seg (60), seg_dur (10 s), fs_lfp (1000 Hz), fpass ([1 200] Hz),
%     tapers ([3 5] = time-bandwidth product, number of tapers)

defaults = struct('n_seg', 60, 'seg_dur', 10, 'fs_lfp', 1000, 'fpass', [1 200], 'tapers', [3 5]);
if ~exist('opts','var'); opts = struct(); end
for name = fieldnames(defaults)'
    if ~isfield(opts, name{1}); opts.(name{1}) = defaults.(name{1}); end
end

chronux_dir = 'C:\Users\Hannah\Dropbox\Matlab\dontaddtopath\chronux_2_11';
addpath(genpath(fullfile(chronux_dir, 'spectral_analysis'))) % only the spectral code: other Chronux folders shadow built-ins

params = struct('Fs', opts.fs_lfp, 'tapers', opts.tapers, 'fpass', opts.fpass, 'pad', 0, 'trialave', 0);
f_info = dir(fbinary);
n_samp = f_info.bytes / (2*nCh);
n_read = round(opts.seg_dur * fs);
starts = round(linspace(fs*10, n_samp - n_read - fs*10, opts.n_seg)); % 0-indexed, skip 10 s at each end

fid = fopen(fbinary, 'r');
cleanup = onCleanup(@() fclose(fid));
S = 0;
for ii = 1:opts.n_seg
    fseek(fid, 2*nCh*starts(ii), 'bof');
    x = double(fread(fid, [nCh n_read], '*int16'));
    x = resample(x(chan_idx,:)' * bit_volts, opts.fs_lfp, fs); % anti-aliased downsampling, samples x channels
    x = x - mean(x);
    [S_seg, f] = mtspectrumc(x, params);
    S = S + S_seg / opts.n_seg;
end
end
