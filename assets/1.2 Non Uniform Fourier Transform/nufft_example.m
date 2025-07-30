%% FINUFFT_example
clear; clc; close all;

%% parameters
fs = 1000;              % sampling frequency
L = 1500;               % length of signal
t_uniform = (0:L-1)/fs; % time vector

%% calculation
% normalize time to [-pi, pi]
T_total = max(t_uniform) - min(t_uniform); 
t_norm = (t_uniform - min(t_uniform)) / T_total * 2*pi - pi;

% create signal
x = 0.7*sin(2*pi*60*t_uniform) + sin(2*pi*120*t_uniform);

% finufft1d1
ms = floor(L/2)+1;
Y = finufft1d1(t_norm, x, -1, 1e-6, ms); 

% frequency domain
tmp = abs(Y/L);
X = tmp(1:floor(L/2)+1);
X(2:end-1) = 2*X(2:end-1);

% frequency axis
f_target = (0:ms-1) * fs / L - fs/4;

%% FFT
Y_fft = fft(x, L);       
P2 = abs(Y_fft/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f_fft = (0:floor(L/2)) * fs / L;

%% create figure
figure('Position', [100, 100, 1200, 600]);

% FINUFFT
subplot(1,2,1);
stem(f_target, X, 'LineWidth', 1.2, 'Marker', 'none');
title('FINUFFT');
xlabel('frequency');
ylabel('|X(f)|');
xlim([0, 200]);       
grid on;
set(gca, 'FontSize', 12);

% FFT
subplot(1,2,2);
stem(f_fft, P1, 'LineWidth', 1.2, 'Marker', 'none', 'Color', [0.8 0.2 0.2]);
title('FFT');
xlabel('frequency');
ylabel('|X(f)|');
xlim([0, 200]);
grid on;
set(gca, 'FontSize', 12);

%% add title
sgtitle('FINUFFT vs FFT', 'FontSize', 16, 'FontWeight', 'bold');