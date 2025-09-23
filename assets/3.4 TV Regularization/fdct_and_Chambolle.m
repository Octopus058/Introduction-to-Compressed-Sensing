clear; clc; close all;

function u = rof_tv_denoise(f, mask, lambda, niter)
    % Chambolle Algorithm
    [N, M] = size(f);
    px = zeros(N, M);  % p
    py = zeros(N, M);
    tau = 0.25;

    for k = 1:niter
        div_p = [px(:,end)-px(:,1), -diff(px,1,2)] + ...
                [py(end,:)-py(1,:); -diff(py,1,1)];

        u = f - lambda * div_p;
        u(mask) = f(mask);

        ux = [diff(u,1,2), u(:,1)-u(:,end)];
        uy = [diff(u,1,1); u(1,:)-u(end,:)];

        px = px + tau * ux;
        py = py + tau * uy;
        norm_p = max(1, sqrt(px.^2 + py.^2));
        px = px ./ norm_p;
        py = py ./ norm_p;
    end
end

% Read the picture
img = im2double(imread('Lena.jpg'));
if size(img,3) == 3
    img = rgb2gray(img);
end
[M, N] = size(img);

% Parameters
sampling_ratio = 0.3;
rng(42);
mask = rand(M, N) < sampling_ratio;

noise_level = 0.01;
b = img .* mask + noise_level * randn(M, N) .* mask;

% Init
x      = b;
maxIter = 200;
lambda  = 0.005;
prev_psnr = 0;
tol_psnr  = 0.002;
lambda_tv = 0.05;
n_tv_iter = 200;

for iter = 1:maxIter
    x_old = x;
    % FDCT
    C = fdct_wrapping(x, 1, 2, ceil(log2(min(M, N)) - 2), 64);
    sigma_n2 = noise_level^2;
    for j = 1:length(C)
        for l = 1:length(C{j})
            coef      = C{j}{l};
            var_coef  = mean(coef(:).^2);
            sigma_x2  = max(var_coef - sigma_n2, eps);
            T         = sigma_n2 / sqrt(sigma_x2);
            C{j}{l}   = sign(coef) .* max(abs(coef) - T, 0);
        end
    end

    % IDFCT and TV regularization
    x_rec = ifdct_wrapping(C, 1, M, N);
    x_rec = rof_tv_denoise(x_rec, mask, lambda_tv, n_tv_iter);

    x = x_rec;
    x(mask) = b(mask);

    % PSNR
    mse = mean((x(:) - img(:)).^2);
    curr_psnr = 10 * log10(1 / mse);
    if mod(iter,10)==0
        fprintf('Iter %d: PSNR = %.2f dB\n', iter, curr_psnr);
    end
    if (iter>1 && abs(curr_psnr - prev_psnr)<tol_psnr) || (curr_psnr-prev_psnr)<0
        fprintf('Stopped at Iter %d\n', iter);
        break;
    end
    prev_psnr = curr_psnr;
end

% Plot the result
figure;
subplot(1,3,1), imshow(img),     title('Origin');
subplot(1,3,2), imshow(b),       title('Sampled');
subplot(1,3,3), imshow(x),       title('Reconstructed');
xlabel(sprintf('PSNR = %.2f dB', curr_psnr), 'FontSize', 10);