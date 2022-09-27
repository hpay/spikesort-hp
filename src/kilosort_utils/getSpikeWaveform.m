function [mean_wave_filt, twave] = getSpikeWaveform(memMap, tspike, Fs,offset_t_final, duration_t_final)
% Inputs
% memMap: raw file
% tspike: spike times
% Fs: samplerate
% offset_t_final: (s)
% duration_t_final: (s)
%
% Outputs
% mean_wave_filt: average filtered waveform 
% twave: timestamps for mean_wave_filt

nCh = memMap.Format{2}(1);

% High pass filter for spikes
Nhp   = 100;  % Order
Fstop = 10;   % Stopband Frequency
Fpass = 800;  % Passband Frequency
b_highpass = firpm(Nhp, [0 Fstop Fpass Fs/2]/(Fs/2), [0 0 1 1]);

% Low pass filter for spikes
Nlp   = 20;   % Order
Fpass = 5000; % Passband Frequency
Fstop = 10000;% Stopband Frequency
b_lowpass = firpm(Nlp, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0]);

% Get the spike waveform
duration_ind_final = round(Fs*duration_t_final);
offset_ind_final = round(Fs*offset_t_final);
duration_ind_initial =  duration_ind_final*2 + (length(b_highpass)+length(b_lowpass))*3; % make sure spike waveform is not affected by the filter
offset_ind_initial = round(duration_ind_initial/2);
start_ind = round(Fs*tspike);
start_ind(start_ind<duration_ind_initial) = [];

% Load the raw data and take the average waveform
% rawdata = readIntanPerChannel(fullfile(raw_folder, raw_file), duration_ind_initial, start_ind-offset_ind_initial);
dataInds = ((start_ind-offset_ind_initial) + (1:duration_ind_initial))';
waveForms = reshape(memMap.Data.x(1:nCh, dataInds(:)), [nCh, duration_ind_initial, length(start_ind)]);
mean_wave = mean(waveForms,3)';

% Filter the average waveform
mean_wave_filt = filtfilt(b_lowpass,1, filtfilt(b_highpass,1, mean_wave));

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
