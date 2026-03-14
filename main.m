%% ========================================================================
%% MULTI-CHANNEL AM BROADCAST & SUPERHETERODYNE RECEIVER SIMULATION
%% ========================================================================
% Description: This script simulates an FDM transmission system and a 
% high-selectivity receiver with Intermediate Frequency (IF) processing.

clear; clc; close all;

%% --- 1. CONFIGURATION & PARAMETERS ---
RadioSystem.Fs_Target = 1e6;       % Global Sampling Rate (1 MHz)
RadioSystem.Freq_IF   = 15000;     % Intermediate Freq (15 kHz)
RadioSystem.Bandwidth = 10e3;     % Signal Bandwidth (10 kHz)
RadioSystem.Carrier1  = 100e3;     % Station A Center Freq
RadioSystem.Carrier2  = 130e3;     % Station B Center Freq

%% --- 2. SOURCE DATA ACQUISITION & ALIGNMENT ---
% Import audio streams
[raw_v1, sample_rate1] = audioread('Short_BBCArabic2.wav');
[raw_v2, sample_rate2] = audioread('Short_FM9090.wav');

% Convert to Mono (Standardize channels)
src1_mono = sum(raw_v1, 2);
src2_mono = sum(raw_v2, 2);

% Calculate Resampling Ratios
up_mod1 = ceil(RadioSystem.Fs_Target / sample_rate1);
up_mod2 = ceil(RadioSystem.Fs_Target / sample_rate2);

% Perform Interpolation to match high-frequency carrier requirements
msg_a_up = interp(src1_mono, up_mod1);
msg_b_up = interp(src2_mono, up_mod2);

% Sync timeline and vector lengths
Fs_Sim = sample_rate1 * up_mod1; 
total_samples = min(length(msg_a_up), length(msg_b_up));
time_vec = (0:total_samples-1) / Fs_Sim;

sig_a = msg_a_up(1:total_samples)';
sig_b = msg_b_up(1:total_samples)';

%% --- 3. TRANSMITTER (FDM MODULATION) ---
% Generate DSB-SC components
mod_ch1 = sig_a .* cos(2 * pi * RadioSystem.Carrier1 * time_vec);
mod_ch2 = sig_b .* cos(2 * pi * RadioSystem.Carrier2 * time_vec);

% Composite Transmitted Signal (Multiplexed)
tx_composite = mod_ch1 + mod_ch2;

% Visualize Frequency Domain of the Composite Signal

points_fft = length(tx_composite);
hz_axis = (-points_fft/2 : points_fft/2-1) * (Fs_Sim/points_fft);
spectrum_tx = abs(fftshift(fft(tx_composite)));

figure('Color', 'w');
plot(hz_axis/1e3, spectrum_tx, 'LineWidth', 1.2);
title('Spectral Density of the Multiplexed Transmission');
xlabel('kHz'); ylabel('Amplitude'); grid on;

%% --- 4. RECEIVER (SUPER-HETERODYNE) ---
% Local Oscillator (High-side injection)
osc_local_freq = RadioSystem.Carrier1 + RadioSystem.Freq_IF;
osc_signal = cos(2 * pi * osc_local_freq * time_vec);

% A. Front-End RF Filtering (Station Selection)
[z_rf, p_rf] = butter(4, [RadioSystem.Carrier1 - RadioSystem.Bandwidth, ...
                          RadioSystem.Carrier1 + RadioSystem.Bandwidth]/(Fs_Sim/2), 'bandpass');
rf_filtered = filter(z_rf, p_rf, tx_composite);

% B. Down-Conversion (Mixing)
mixed_signal = rf_filtered .* osc_signal;

% C. Intermediate Frequency (IF) Stage
[z_if, p_if] = butter(4, [RadioSystem.Freq_IF - RadioSystem.Bandwidth, ...
                          RadioSystem.Freq_IF + RadioSystem.Bandwidth]/(Fs_Sim/2), 'bandpass');
if_filtered = filter(z_if, p_if, mixed_signal);

% D. Synchronous Detector (Baseband Recovery)
detector_osc = cos(2 * pi * RadioSystem.Freq_IF * time_vec);
baseband_raw = if_filtered .* detector_osc;

[z_audio, p_audio] = butter(4, (RadioSystem.Bandwidth)/(Fs_Sim/2), 'low');
recovered_audio = filter(z_audio, p_audio, baseband_raw);

%% --- 5. OUTPUT GENERATION ---
% Decimate back to original audio sampling rate
final_stream = downsample(recovered_audio, up_mod1);

% Export Processed Audio
audiowrite('Recovered_Station1.wav', final_stream, sample_rate1);
disp('>>> Receiver processing complete. File saved.');

%% --- 6. ADVANCED SYSTEM ANALYSIS ---
% Helper for FFT Calculation
calc_fft = @(data) abs(fftshift(fft(data)));

figure('Name', 'System Stage Diagnostics');
subplot(3,1,1);
plot(hz_axis/1e3, calc_fft(rf_filtered));
title('Stage 1: RF Filter Output (Isolated Channel)');
xlim([50 150]); grid on;

subplot(3,1,2);
plot(hz_axis/1e3, calc_fft(if_filtered));
title('Stage 2: IF Filter Output (Down-converted)');
xlim([-30 30]); grid on;

subplot(3,1,3);
plot(hz_axis/1e3, calc_fft(baseband_raw));
title('Stage 3: Post-Detector Spectrum (Prior to Audio LPF)');
xlim([-40 40]); grid on;

%% --- 7. STRESS TESTING (OFFSETS & INTERFERENCE) ---
% Scenario 1: LO Frequency Drift (Task 5)
drift_values = [100, 1000]; % 0.1kHz and 1kHz
for drift = drift_values
    drift_baseband = if_filtered .* cos(2 * pi * (RadioSystem.Freq_IF + drift) * time_vec);
    drift_audio = filter(z_audio, p_audio, drift_baseband);
    
    out_drift = downsample(drift_audio, up_mod1);
    fname = sprintf('Drift_Simulation_%dHz.wav', drift);
    audiowrite(fname, out_drift / max(abs(out_drift)), sample_rate1);
    disp('>>> Testing Offsets Audio processing complete. File saved.');
end

% Scenario 2: Image Channel Interference (Task 4 - Disable RF Stage)
leaked_mix = tx_composite .* osc_signal; 
leaked_if = filter(z_if, p_if, leaked_mix);
leaked_baseband = leaked_if .* detector_osc;
leaked_audio = filter(z_audio, p_audio, leaked_baseband);

out_corrupted = downsample(leaked_audio, up_mod1);
audiowrite('Image_Interference_Test.wav', out_corrupted / max(abs(out_corrupted)), sample_rate1);
    disp('>>> Image_Interference_Test audio processing complete. File saved.');