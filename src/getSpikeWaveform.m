function [mean_wave_filt, twave] = getSpikeWaveform(memMap, tspike, Fs,offset_t_final, duration_t_final, hp_pass)
% Inputs
% memMap: raw file
% tspike: spike times
% Fs: samplerate
% offset_t_final: (s)
% duration_t_final: (s)
% hp_pass (optional): high-pass cutoff (Hz), default 300
%
% Outputs
% mean_wave_filt: average filtered waveform
% twave: timestamps for mean_wave_filt
%
% Filtering: 3rd-order Butterworth high-pass at 300 Hz (Kilosort 4's filter) and 3rd-order Butterworth
% low-pass at 5 kHz, both forward and backward (zero phase), with 5 ms of padding on each side. Replaced a 100-tap equiripple FIR high-pass
% (passband from 800 Hz) whose end taps added a small bump ~1.6 ms from the spike.

nCh = memMap.Format{2}(1);

hp_cut = 300; if exist('hp_pass', 'var') && ~isempty(hp_pass); hp_cut = hp_pass; end
lp_cut = 5000;
[b_hp, a_hp] = butter(3, hp_cut/(Fs/2), 'high');
[b_lp, a_lp] = butter(3, lp_cut/(Fs/2), 'low');

% Get the spike waveform
duration_ind_final = round(Fs*duration_t_final);
offset_ind_final = round(Fs*offset_t_final);
duration_ind_initial =  duration_ind_final*2 + 2*round(0.005*Fs); % 5 ms extra on each side so the filters settle before the cropped window (changes the kept window by <0.01%)
offset_ind_initial = round(duration_ind_initial/2);
start_ind = round(Fs*tspike);
start_ind(start_ind<duration_ind_initial) = [];

% Load the raw data and take the average waveform
% fread each window directly (~10x faster than indexing the memmapfile; same result)
fid = fopen(memMap.Filename, 'r');
wave_sum = zeros(nCh, duration_ind_initial);
for first_samp = reshape(start_ind - offset_ind_initial, 1, []) % 0-indexed first sample of each window
    fseek(fid, 2*nCh*first_samp, 'bof');
    wave_sum = wave_sum + double(fread(fid, [nCh duration_ind_initial], '*int16'));
end
fclose(fid);
mean_wave = (wave_sum / numel(start_ind))';

% Filter the average waveform
mean_wave_filt = filtfilt(b_lp, a_lp, filtfilt(b_hp, a_hp, mean_wave));

% % Upsample (for more precise estimate of duration?)
% us = 10;
% Fs_upsample = Fs*10;
% mean_wave_upsample = interpft(mean_wave_filt,length(mean_wave_filt)*10);
% mean_wave_upsample = mean_wave_upsample((offset_ind_initial-offset_ind_final)*us + (1:duration_ind_final*us),:);
% twave_upsample = (1:length(mean_wave_upsample))/(Fs*us)-offset_t_final;

% Crop to remove edges affected by filter
wave_mask = (offset_ind_initial-offset_ind_final) + (1:duration_ind_final);
mean_wave_filt = mean_wave_filt(wave_mask,:);

% Create timestamps for plotting purposes
twave = (1:nnz(wave_mask))/Fs-offset_t_final;
