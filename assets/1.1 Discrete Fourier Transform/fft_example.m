%% FFT_example
clear; clc; close all;

%% parameters
Fs = 1000;            % sampling frequency
T = 1/Fs;             % sampling period
L = 1500;             % length of signal
t = (0:L-1)*T;        % time vector

S = 0.7 * sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));

%% time domain
figure('Name','FFT','NumberTitle','off')
subplot(2,1,1)
plot(t,X,'b','LineWidth',1.2)
title('time domain')
xlabel('time')
ylabel('amplitude')
grid on
xlim([0 0.2])

%% FFT
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% frequency axis
f = Fs*(0:(L/2))/L;

%% frequency domain
subplot(2,1,2)
plot(f,P1,'r','LineWidth',1.5)
title('frequency domain')
xlabel('frequency')
ylabel('amplitude')
grid on
xlim([0 250])
