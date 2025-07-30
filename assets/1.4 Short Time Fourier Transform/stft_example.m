%% STFT_example
clear; clc; close all;

%% 参数
Fs = 1000;            % sampling frequency
T = 1/Fs;             % sampling period
L = 2000;             % length of signal
t = (0:L-1)*T;        % time vector

% create signal
X1 = [sin(2*pi*20*t(1:L/2)), 2*sin(2*pi*10*t(L/2+1:end))];
X2 = [2*sin(2*pi*10*t(1:L/2)), sin(2*pi*20*t(L/2+1:end))];

%% create figure
figure('Name','STFT','NumberTitle','off', 'Position', [100, 100, 1200, 600])

% window function
window_size = 0.5; 
ts = [0.5, 1.0, 1.5, 2.0]; % t_s
y_positions = [0, 0.5, 1.0, 1.5];
colors = {'b', 'g', 'r', 'm'};
freq_range = [0, 50];

%% STFT of x1
subplot(1,2,1)
hold on; grid on;
title('X1')
xlabel('frequency')
ylabel('ts')
zlabel('amplitude')
view(45, 30)

for i = 1:length(ts)
    t_end = ts(i);
    t_start = t_end - window_size;
    
    idx_start = find(t >= t_start, 1);
    idx_end = find(t <= t_end, 1, 'last');
    segment = X1(idx_start:idx_end);
    seg_length = length(segment);
    
    % FFT
    Y = fft(segment);
    P2 = abs(Y/seg_length);
    P1 = P2(1:floor(seg_length/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % frequency axis
    f = Fs*(0:(floor(seg_length/2)))/seg_length;
    
    freq_idx = find(f >= freq_range(1) & f <= freq_range(2));
    f_display = f(freq_idx);
    P1_display = P1(freq_idx);
    
    plot3(f_display, repmat(y_positions(i), size(f_display)), P1_display, ...
          'Color', colors{i}, 'LineWidth', 2)
end

% axis
xlim(freq_range)
ylim([-0.1 1.6])
zlim([0 2.5])
set(gca, 'FontSize', 12)
legend('0s', '0.5s', '1s', '1.5s', 'Location', 'northeast')

%% STFT of x2
subplot(1,2,2)
hold on; grid on;
title('X2')
xlabel('frequency')
ylabel('ts')
zlabel('amplitude')
view(45, 30)

for i = 1:length(ts)
    t_end = ts(i);
    t_start = t_end - window_size;
    
    idx_start = find(t >= t_start, 1);
    idx_end = find(t <= t_end, 1, 'last');
    segment = X2(idx_start:idx_end);
    seg_length = length(segment);
    
    % FFT
    Y = fft(segment);
    P2 = abs(Y/seg_length);
    P1 = P2(1:floor(seg_length/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % frequency axis
    f = Fs*(0:(floor(seg_length/2)))/seg_length;
    
    freq_idx = find(f >= freq_range(1) & f <= freq_range(2));
    f_display = f(freq_idx);
    P1_display = P1(freq_idx);
    
    plot3(f_display, repmat(y_positions(i), size(f_display)), P1_display, ...
          'Color', colors{i}, 'LineWidth', 2)
end

% axis
xlim(freq_range)
ylim([-0.1 1.6])
zlim([0 2.5])
set(gca, 'FontSize', 12)
legend('0s', '0.5s', '1s', '1.5s', 'Location', 'northeast')

% title
sgtitle('STFT', 'FontSize', 16, 'FontWeight', 'bold')