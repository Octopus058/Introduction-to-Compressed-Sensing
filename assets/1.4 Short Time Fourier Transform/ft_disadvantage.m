%% FT_disadvantage
clear; clc; close all;

%% parameters
Fs = 1000;            % sampling frequency
T = 1/Fs;             % sampling period
L = 2000;             % length of signal
t = (0:L-1)*T;        % time vector

% create signal
X1 = [sin(2*pi*20*t(1:L/2)), 2*sin(2*pi*10*t(L/2+1:end))];
X2 = [2*sin(2*pi*10*t(1:L/2)), sin(2*pi*20*t(L/2+1:end))];

%% create figure
figure('Name','FT`s disadvantage','NumberTitle','off', 'Position', [100, 100, 1200, 600])

%% time domain
subplot(1,2,1)
plot(t, X1, 'b', 'LineWidth', 1.5)
hold on
plot(t, X2, 'r-', 'LineWidth', 1.5)
title('time domain')
xlabel('time')
ylabel('amplitude')
grid on
xlim([0 2])
ylim([-2.5 2.5])

% 1s
line([1 1], [-2.5 2.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.2)
text(0.5, 2.3, 'Previous second', 'HorizontalAlignment', 'center', 'FontSize', 12, 'BackgroundColor', 'w')
text(1.5, 2.3, 'Next second', 'HorizontalAlignment', 'center', 'FontSize', 12, 'BackgroundColor', 'w')

% legend
legend('X1: Previous second(20Hz), Next second(10Hz)', ...
       'X2: Previous second(10Hz), Next second(20Hz)', ...
       'Location', 'southwest')
set(gca, 'FontSize', 12)

%% frequency domain
subplot(1,2,2)

% FFT of x1
Y1 = fft(X1);
P2_X1 = abs(Y1/L);
P1_X1 = P2_X1(1:L/2+1);
P1_X1(2:end-1) = 2*P1_X1(2:end-1);

% FFT of x2
Y2 = fft(X2);
P2_X2 = abs(Y2/L);
P1_X2 = P2_X2(1:L/2+1);
P1_X2(2:end-1) = 2*P1_X2(2:end-1);

% frequency axis
f = Fs*(0:(L/2))/L;

% frequency domain
plot(f, P1_X1, 'b', 'LineWidth', 1.5)
hold on
plot(f, P1_X2, 'r--', 'LineWidth', 1.5)
title('frequency domain')
xlabel('frequency')
ylabel('amplitude')
grid on
xlim([0 50])
set(gca, 'FontSize', 12)

% legend
legend('X1', 'X2', 'Location', 'northeast')

% title
sgtitle('FT`s disadvantage', 'FontSize', 16, 'FontWeight', 'bold')