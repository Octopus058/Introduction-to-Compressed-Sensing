%% CoSaMP_example
clear; clc; close all;
rng(42);

%% Parameters
N = 4000;              % signal length
Fs = 1000;             % sampling frequency
T = N/Fs;              % signal duration
t = (0:N-1)/Fs;        % time vector

sparsity_level = 8;    % target sparsity
loss_percentage = 0.97;% data loss ratio
noise_level = 0.1;     % noise level

%% Signal Generation
freqs = [10, 25, 40, 60];
amps  = [1.0, 0.8, 0.6, 0.4];
x_original = zeros(1, N);
for i = 1:length(freqs)
    x_original = x_original + amps(i)*sin(2*pi*freqs(i)*t);
    if mod(i,2)==0
        x_original = x_original + 0.7*amps(i)*cos(2*pi*freqs(i)*t);
    end
end
x_original = x_original / max(abs(x_original));
x_original = x_original(:);

%% Data Loss and Noise
mask       = rand(N,1) > loss_percentage;
x_observed = x_original .* mask + noise_level*randn(N,1);

%% Construct Observation System
obs_idx = find(mask);
M       = length(obs_idx);
y       = x_observed(obs_idx);
A = zeros(M, N);
for i = 1:M
    A(i, obs_idx(i)) = 1;
end

%% Construct DFT Sparse Basis and Sensing Matrix
Psi = dftmtx(N) / sqrt(N);
for i = 1:N
    Psi(:,i) = Psi(:,i) / norm(Psi(:,i));
end
Phi = A * Psi;

%% CoSaMP
max_iter    = 10;
K           = sparsity_level;
residual    = y;
theta_hat   = zeros(N,1);
support_set = [];
res_history = zeros(max_iter,1);

for iter = 1:max_iter
    % 1. Residual correlations
    correlations = abs(Phi' * residual);
    % 2. Select 2K atoms
    [~, idx]      = sort(correlations,'descend');
    candidate_set = idx(1:min(2*K, length(idx)));
    % 3. Merge support set
    merged_set    = union(support_set, candidate_set);
    % 4. Least squares estimation
    Phi_sub = Phi(:, merged_set);
    x_ls    = Phi_sub \ y;
    % 5. Keep K largest coefficients
    [~, sidx]   = sort(abs(x_ls),'descend');
    support_set = merged_set(sidx(1:min(K, length(x_ls))));
    % 6. Update sparse coefficients
    theta_hat(support_set) = Phi(:, support_set) \ y;
    % 7. Update residual
    residual = y - Phi(:, support_set) * theta_hat(support_set);
    res_history(iter) = norm(residual);
    if res_history(iter) < 1e-6
        break;
    end
end

%% Reconstruct Signal
x_reconstructed = Psi * theta_hat;
x_final         = real(x_reconstructed);
x_final = medfilt1(x_final, 7);

%% Result Visualization
figure('Position',[100,100,900,800],'Name','CoSaMP Reconstruction Results');

% Time domain
subplot(4,1,1);
plot(t, x_original,'b','LineWidth',1.8); hold on;
stem(t, x_observed, 'r','Marker','o','MarkerSize',3,'MarkerFaceColor','r');
title(['Original Signal and Observations (Data loss: ',num2str(loss_percentage*100),'%, Noise level: ',num2str(noise_level),')']);
legend('Original signal','Observed samples'); xlabel('Time'); ylabel('Amplitude');
xlim([0,T]); grid on;

subplot(4,1,2);
plot(t, x_original,'b','LineWidth',1.8); hold on;
plot(t, x_final,'r--','LineWidth',1.5);
title('Original vs Reconstructed (Time domain)');
legend('Original signal','CoSaMP'); xlabel('Time'); ylabel('Amplitude');
xlim([0,T]); grid on;

% Reconstruction error
subplot(4,1,3);
error_signal = x_original - x_final;
plot(t, error_signal,'g','LineWidth',1.5);
title('Reconstruction Error'); xlabel('Time'); ylabel('Error');
xlim([0,T]); grid on;

% Frequency domain
subplot(4,1,4);
halfN = floor(N/2) + 1;
freq  = linspace(0, Fs/2, halfN);
Y_o   = fft(x_original);
P_o   = abs(Y_o(1:halfN)/N);
P_o(2:end-1) = 2*P_o(2:end-1);
Y_r   = fft(x_final);
P_r   = abs(Y_r(1:halfN)/N);
P_r(2:end-1) = 2*P_r(2:end-1);
P_o_db = 20*log10(P_o + eps);
P_r_db = 20*log10(P_r + eps);
plot(freq, P_o_db,'b','LineWidth',1.5); hold on;
plot(freq, P_r_db,'r--','LineWidth',1.2);
title('Original vs Reconstructed (Frequency domain)');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Original signal','Reconstructed'); xlim([0,100]); grid on;
ylim([-50, 0]);

%% Performance Metrics
error_final     = x_original - x_final;
rel_error_final = norm(error_final)/norm(x_original);

fprintf('===== CoSaMP Reconstruction Result Statistics =====\n');
fprintf('Signal length: %d\n', N);
fprintf('Sampling frequency: %d Hz\n', Fs);
fprintf('Data loss percentage: %.1f%%\n', loss_percentage*100);
fprintf('Number of observations: %d (%.1f%%)\n', M, (M/N)*100);
fprintf('Sparse basis: DFT\n');
fprintf('Target sparsity: %d\n', sparsity_level);
fprintf('CoSaMP iterations: %d\n', iter);
fprintf('Relative reconstruction error: %.4f\n', rel_error_final);