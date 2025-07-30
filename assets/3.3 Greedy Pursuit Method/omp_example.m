%% OMP_example
clear; clc; close all;
rng(42);

%% Parameters
N = 500;               % signal length
Fs = 1000;             % sampling frequency
T = N/Fs;              % signal duration
t = (0:N-1)/Fs;        % time vector

sparsity_level = 8;    % estimated sparsity
loss_percentage = 0.3; % 30% data loss
noise_level = 0.1;     % noise level

%% Signal Generation
% frequency components
freqs = [10, 25, 40, 60];
amps = [1.0, 0.8, 0.6, 0.4];

x_original = zeros(1, N);
for i = 1:length(freqs)
    x_original = x_original + amps(i) * sin(2*pi*freqs(i)*t);
    if mod(i, 2) == 0
        x_original = x_original + 0.7*amps(i) * cos(2*pi*freqs(i)*t);
    end
end

% normalization
x_original = x_original / max(abs(x_original));
x_original = x_original(:);

%% Data Loss and Noise
mask = rand(N, 1) > loss_percentage;
x_observed = x_original .* mask + noise_level*randn(N,1);

%% Construct Observation System
obs_idx = find(mask);
M       = length(obs_idx);
y       = x_observed(obs_idx);
A = zeros(M, N);
for i = 1:M
    A(i, obs_idx(i)) = 1;
end

%% Construct DCT Sparse Basis and Sensing Matrix
Psi = dctmtx(N)';
% normalize basis vectors
for i = 1:N
    Psi(:, i) = Psi(:, i) / norm(Psi(:, i));
end

% sensing matrix
Phi = A * Psi;
% normalize columns
col_norms = sqrt(sum(Phi.^2, 1));
col_norms(col_norms < eps) = 1;
Phi = Phi ./ col_norms;

%% OMP
max_iter = 33;          % maximum iterations
residual = y;           % initial residual
index_set = [];         % support set
theta_hat = zeros(size(Psi, 2), 1);  % sparse coefficients estimate
res_history = zeros(max_iter, 1);    % residual history
selected_atoms = false(size(Psi, 2), 1); % track selected atoms

% OMP main loop
for iter = 1:max_iter
    % compute projections of residual onto sensing matrix
    projections = abs(Phi' * residual);
    projections(selected_atoms) = -inf; % exclude selected atoms
    % find index of max projection
    [max_val, new_idx] = max(projections);
    % check if atom found
    if max_val < 1e-6 || isinf(max_val)
        break;
    end
    % mark and add new atom
    selected_atoms(new_idx) = true;
    index_set = [index_set, new_idx];
    % build submatrix with selected atoms
    Phi_sub = Phi(:, selected_atoms);
    % least squares solution
    if iter == 1
        Q = Phi_sub / norm(Phi_sub);
        R = 1;
        x_ls = Q' * y;
    else
        [Q, R] = qr(Phi_sub, 0);
        x_ls = R \ (Q' * y);
    end
    % update residual
    residual = y - Phi_sub * x_ls;
    res_history(iter) = norm(residual);
    % check stopping criteria
    if norm(residual) < 1e-4 || iter >= max_iter
        break;
    end
end

% build sparse coefficient estimate
theta_hat(selected_atoms) = x_ls;

%% Reconstruct Signal
x_reconstructed = Psi * theta_hat;
x_final = x_reconstructed;
x_final(obs_idx) = x_observed(obs_idx);

%% Result Visualization
figure('Position', [100, 100, 900, 800], 'Name', 'OMP Reconstruction Results');

% Time domain
subplot(4,1,1);
plot(t, x_original, 'b', 'LineWidth', 1.8);
hold on;
stem(t, x_observed, 'r', 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
title(['Original and Observed Signals (Data loss: ', num2str(loss_percentage*100), '%, Noise level: ', num2str(noise_level), ')']);
legend('Original signal', 'Observed samples', 'Location', 'best');
xlabel('Time');
ylabel('Amplitude');
xlim([0, T]);
grid on;

subplot(4,1,2);
plot(t, x_original, 'b', 'LineWidth', 1.8);
hold on;
plot(t, x_final, 'r--', 'LineWidth', 1.5);
title('Original vs Reconstructed (Time domain)');
legend('Original signal', 'OMP', 'Location', 'best');
xlabel('Time');
ylabel('Amplitude');
xlim([0, T]);
grid on;

% Reconstruction error
subplot(4,1,3);
error_signal = x_original - x_final;
plot(t, error_signal, 'g', 'LineWidth', 1.5);
title('Reconstruction Error');
xlabel('Time');
ylabel('Error');
xlim([0, T]);
grid on;

% Frequency domain
subplot(4,1,4);
freq = linspace(0, Fs/2, floor(N/2)+1);
Y_orig = fft(x_original);
P_orig = abs(Y_orig(1:floor(N/2)+1)/N);
Y_rec = fft(x_final);
P_rec = abs(Y_rec(1:floor(N/2)+1)/N);

plot(freq, 20*log10(P_orig), 'b', 'LineWidth', 1.5);
hold on;
plot(freq, 20*log10(P_rec), 'r--', 'LineWidth', 1.2);
title('Original vs Reconstructed (Frequency domain)');
xlabel('Frequency');
ylabel('Magnitude (dB)');
legend('Original signal', 'Reconstructed', 'Location', 'best');
xlim([0, 100]);
grid on;

%% Compute Reconstruction Metrics
% final error
error_final = x_original - x_final;
rel_error_final = norm(error_final) / norm(x_original);

fprintf('===== Reconstruction Result Statistics =====\n');
fprintf('Signal length: %d\n', N);
fprintf('Sampling frequency: %d Hz\n', Fs);
fprintf('Data loss percentage: %.1f%%\n', loss_percentage*100);
fprintf('Number of observations: %d (%.1f%%)\n', M, (M/N)*100);
fprintf('Sparse basis: DCT\n');
fprintf('Estimated sparsity: %d\n', sparsity_level);
fprintf('OMP iterations: %d\n', iter);
fprintf('Relative reconstruction error: %.4f\n', rel_error_final);