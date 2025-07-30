%% DWT_example
clear; clc; close all;

%% Parameters
fs = 1000;              % sampling frequency
N = 64;                 % number of samples
t = (0:N-1)/fs;         % time vector

% Original signal
f1 = 50;
f2 = 100;
x_clean = 2*sin(2*pi*f1*t) + sin(2*pi*f2*t);

% Gaussian white noise
SNR = 10;
noise_power = var(x_clean) / (10^(SNR/10));
noise = sqrt(noise_power) * randn(1, N);
x_noisy = x_clean + noise;

%% Wavelet decomposition
wavelet = 'db4';        % db4 wavelet
level = 4;              % 4 levels of decomposition

[c, l] = wavedec(x_noisy, level, wavelet);

% Reconstruct approximation and detail coefficients at each level
A = cell(1, level);
D = cell(1, level);
for k = 1:level
    A{k} = wrcoef('a', c, l, wavelet, k);  % Approximation at level k
    D{k} = wrcoef('d', c, l, wavelet, k);  % Detail at level k
end

% Reconstruction
x_recon = waverec(c, l, wavelet);

% Compute frequency bands for each sub-band
freq_bands = cell(1, level+1);
for k = 1:level
    high_freq = fs / (2^k);
    low_freq  = fs / (2^(k+1));
    freq_bands{level-k+2} = sprintf('D%d: %.1f–%.1f Hz', k, low_freq, high_freq);
end
freq_bands{1} = sprintf('A%d: 0–%.1f Hz', level, fs/(2^(level+1)));

%% create figure
figure('Position', [100, 100, 800, 900]);
numPlots = level + 2;

% 1. Noisy signal
subplot(numPlots, 1, 1);
plot(t, x_noisy, 'color', [0.6 0 0.6], 'LineWidth', 1.5);
title(sprintf('Noisy Signal (SNR=%ddB, N=%d)', SNR, N));
xlabel('time'); ylabel('amplitude');
xlim([t(1), t(end)]); grid on; ylim([-5, 5]);

% 2. Approximation factor
subplot(numPlots, 1, 2);
plot(t, A{level}, 'r', 'LineWidth', 1.5);
title(freq_bands{1});
xlabel('time'); ylabel('amplitude');
xlim([t(1), t(end)]); grid on; ylim([-3, 3]);

% 3. Detail factor
for k = level:-1:1
    subplot(numPlots, 1, level - k + 2);
    plot(t, D{k}, 'color', [0 0.5 0], 'LineWidth', 1.2);
    title(freq_bands{level - k + 2});
    xlabel('time'); ylabel('amplitude');
    xlim([t(1), t(end)]); grid on;
    % Adjust y-axis based on subband
    if k == 4 || k == 3
        ylim([-3, 3]);
    else
        ylim([-1.5, 1.5]);
    end
end

% 4. Reconstructed signal
subplot(numPlots, 1, numPlots);
plot(t, x_recon, 'b', 'LineWidth', 1.5);
hold on;
plot(t, x_noisy, 'r--', 'LineWidth', 1.2);
xlabel('time'); ylabel('amplitude');
xlim([t(1), t(end)]); grid on; ylim([-5, 5]);
legend('Reconstructed Signal', 'Noisy Signal', 'Location', 'best');
title('Reconstructed Signal');

sgtitle('DWT', 'FontSize', 12);