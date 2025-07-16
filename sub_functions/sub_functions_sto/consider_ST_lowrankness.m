clear;
close all;

addpath(genpath('sub_functions'))


%% Selecting noise condition
deg.gaussian_sigma      = 0.1; % standard derivation of Gaussian noise
deg.sparse_rate         = 0.05;
deg.stripe_rate         = 0.05;
deg.stripe_intensity    = 0.5;

% image = 'JasperRidge';
% image = 'PaviaU120';
% image = 'WashingtonDC';
image = 'Salinas';


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

% HSI_clean = single(HSI_clean);
% HSI_noisy = single(HSI_noisy);


%% Setting difference operator
D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
% Dv      = @(z) z([2:end, 1],:,:) - z;
% Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Ds      = @(z) z(:, :, [2:end, 1], :) - z;
Dst     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);


%% Setting common parameters
rho = 0.95;
boundary = 'circulant';
blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1];
stopcri_index = 5;

name_method = 'S3TTV';
name_params_savetext = ...
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)); % S3TTV


%% Loading restored HS images
save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_tint', num2str(deg.stripe_intensity), '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_restored' ...
);

HSI_restored = double(HSI_restored);

HSI_rand     = rand(hsi.sizeof);

%% Generating structure tensors
DU_true   = D(HSI_clean);
DV        = D(HSI_noisy);
DU        = D(HSI_restored);
DU_rand   = D(HSI_rand);

DsU_true   = Ds(HSI_clean);
DsV        = Ds(HSI_noisy);
DsU        = Ds(HSI_restored);
DsU_rand   = Ds(HSI_rand);

DDsU_true   = D(DsU_true);
DDsV        = D(DsV);
DDsU        = D(DsU);
DDsU_rand   = D(DsU_rand);


ST_U_true   = reshape(HSI_clean, [hsi.n1*hsi.n2, hsi.n3]);
ST_V        = reshape(HSI_noisy, [hsi.n1*hsi.n2, hsi.n3]);
ST_U        = reshape(HSI_restored, [hsi.n1*hsi.n2, hsi.n3]);
ST_U_rand   = reshape(HSI_rand, [hsi.n1*hsi.n2, hsi.n3]);

ST_DU_true   = reshape(DU_true, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DV        = reshape(DV, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DU        = reshape(DU, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DU_rand   = reshape(DU_rand, [hsi.n1*hsi.n2, hsi.n3*2]);

ST_DsU_true   = reshape(DsU_true, [hsi.n1*hsi.n2, hsi.n3]);
ST_DsV        = reshape(DsV, [hsi.n1*hsi.n2, hsi.n3]);
ST_DsU        = reshape(DsU, [hsi.n1*hsi.n2, hsi.n3]);
ST_DsU_rand   = reshape(DsU_rand, [hsi.n1*hsi.n2, hsi.n3]);

ST_DDsU_true   = reshape(DDsU_true, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DDsV        = reshape(DDsV, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DDsU        = reshape(DDsU, [hsi.n1*hsi.n2, hsi.n3*2]);
ST_DDsU_rand   = reshape(DDsU_rand, [hsi.n1*hsi.n2, hsi.n3*2]);


%% Plotting ST_U
% is_plot_ST_U = 1;
if exist('is_plot_ST_U', 'var') && is_plot_ST_U
figure
subplot(1, 3, 1)
imagesc(ST_U_true)
title('Ground-truth')

subplot(1, 3, 2)
imagesc(ST_V)
title('Observation')

subplot(1, 3, 3)
imagesc(ST_U)
title('S_{3}TTV')

sgtitle('ST U');
end


%% Plotting ST_DsU
% is_plot_ST_DsU = 1;
if exist('is_plot_ST_DsU', 'var') && is_plot_ST_DsU
figure
subplot(1, 3, 1)
imagesc(ST_DsU_true)
title('Ground-truth')

subplot(1, 3, 2)
imagesc(ST_DsV)
title('Observation')

subplot(1, 3, 3)
imagesc(ST_DsU)
title('S_{3}TTV')

sgtitle('ST DsU');

end


%% Plotting ST_DU
% is_plot_ST_DU = 1;
if exist('is_plot_ST_DU', 'var') && is_plot_ST_DU
figure
subplot(1, 3, 1)
imagesc(ST_DU_true)
title('Ground-truth')

subplot(1, 3, 2)
imagesc(ST_V)
title('Observation')

subplot(1, 3, 3)
imagesc(ST_DU)
title('S_{3}TTV')

sgtitle('ST DU');

end

%% Plotting ST_DDsU
% is_plot_ST_DDsU = 1;
if exist('is_plot_ST_DDsU', 'var') && is_plot_ST_DDsU
figure
subplot(1, 3, 1)
imagesc(ST_DDsU_true)
title('Ground-truth')

subplot(1, 3, 2)
imagesc(ST_DDsV)
title('Observation')

subplot(1, 3, 3)
imagesc(ST_DDsU)
title('S_{3}TTV')

sgtitle('ST DDsU');

end


%% Plotting singular values for U, DsU, and DDsU
is_plot_SVD_3types = 1;
if exist('is_plot_SVD_3types', 'var') && is_plot_SVD_3types

% 特異値分解
[~, S_ST_U, ~] = svd(ST_U);
[~, S_ST_U_true, ~] = svd(ST_U_true);

[~, S_ST_DsU, ~] = svd(ST_DsU);
[~, S_ST_DsU_true, ~] = svd(ST_DsU_true);

[~, S_ST_DDsU, ~] = svd(ST_DDsU);
[~, S_ST_DDsU_true, ~] = svd(ST_DDsU_true);


% 特異値を取得
singular_values_ST_U = diag(S_ST_U);
singular_values_ST_U_true = diag(S_ST_U_true);

singular_values_ST_DsU = diag(S_ST_DsU);
singular_values_ST_DsU_true = diag(S_ST_DsU_true);

singular_values_ST_DDsU = diag(S_ST_DDsU);
singular_values_ST_DDsU_true = diag(S_ST_DDsU_true);


% 特異値を対数スケールでプロット
% figure(4);
% subplot(1, 3, 1)
% semilogy(singular_values1, 'o-', 'DisplayName', 'S3TTV');
% hold on;
% semilogy(singular_values2, 'x-', 'DisplayName', 'Ground-truth');
% xlabel('Index');
% ylabel('Singular Value (log scale)');
% title('ST U');
% legend;
% grid on;
% hold off;


% 累積寄与率の計算
explained_variance_ST_U = cumsum(singular_values_ST_U) / sum(singular_values_ST_U);
explained_variance_ST_U_true = cumsum(singular_values_ST_U_true) / sum(singular_values_ST_U_true);

explained_variance_ST_DsU = cumsum(singular_values_ST_DsU) / sum(singular_values_ST_DsU);
explained_variance_ST_DsU_true = cumsum(singular_values_ST_DsU_true) / sum(singular_values_ST_DsU_true);

explained_variance_ST_DDsU = cumsum(singular_values_ST_DDsU) / sum(singular_values_ST_DDsU);
explained_variance_ST_DDsU_true = cumsum(singular_values_ST_DDsU_true) / sum(singular_values_ST_DDsU_true);


% 累積寄与率をプロット
figure(4);
subplot(1, 3, 1)
plot(explained_variance_ST_U, 'o-', 'DisplayName', 'S3TTV');
hold on;
plot(explained_variance_ST_U_true, 'x-', 'DisplayName', 'Ground-truth');
xlabel('Number of Singular Values');
ylabel('Cumulative Explained Variance');
title('ST U');
legend;
grid on;
hold off;

subplot(1, 3, 2)
plot(explained_variance_ST_DsU, 'o-', 'DisplayName', 'S3TTV');
hold on;
plot(explained_variance_ST_DsU_true, 'x-', 'DisplayName', 'Ground-truth');
xlabel('Number of Singular Values');
ylabel('Cumulative Explained Variance');
title('ST DsU');
legend;
grid on;
hold off;

subplot(1, 3, 3)
plot(explained_variance_ST_DDsU, 'o-', 'DisplayName', 'S3TTV');
hold on;
plot(explained_variance_ST_DDsU_true, 'x-', 'DisplayName', 'Ground-truth');
xlabel('Number of Singular Values');
ylabel('Cumulative Explained Variance');
title('ST DDsU');
legend;
grid on;
hold off;

sgtitle('explained variance');

end

%% Plotting singular values for DDsU, DDsV, and DDsU_rand
% is_plot_SVD_DDsU = 1;
if exist('is_plot_SVD_DDsU', 'var') && is_plot_SVD_DDsU

% 特異値分解
[~, S_ST_U, ~] = svd(ST_U);
[~, S_ST_U_true, ~] = svd(ST_U_true);
[~, S_ST_U_rand, ~] = svd(ST_U_rand);

[~, S_ST_DDsU, ~] = svd(ST_DDsU);
[~, S_ST_DDsU_true, ~] = svd(ST_DDsU_true);
[~, S_ST_DDsU_rand, ~] = svd(ST_DDsU_rand);


% 特異値を取得
singular_values_ST_U = diag(S_ST_U);
singular_values_ST_U_true = diag(S_ST_U_true);
singular_values_ST_U_rand = diag(S_ST_U_rand);

singular_values_ST_DDsU = diag(S_ST_DDsU);
singular_values_ST_DDsU_true = diag(S_ST_DDsU_true);
singular_values_ST_DDsU_rand = diag(S_ST_DDsU_rand);


% 特異値を対数スケールでプロット
% figure(4);
% subplot(1, 3, 1)
% semilogy(singular_values1, 'o-', 'DisplayName', 'S3TTV');
% hold on;
% semilogy(singular_values2, 'x-', 'DisplayName', 'Ground-truth');
% xlabel('Index');
% ylabel('Singular Value (log scale)');
% title('ST U');
% legend;
% grid on;
% hold off;


% 累積寄与率の計算
explained_variance_ST_U = cumsum(singular_values_ST_U) / sum(singular_values_ST_U);
explained_variance_ST_U_true = cumsum(singular_values_ST_U_true) / sum(singular_values_ST_U_true);
explained_variance_ST_U_rand = cumsum(singular_values_ST_U_rand) / sum(singular_values_ST_U_rand);


explained_variance_ST_DDsU = cumsum(singular_values_ST_DDsU) / sum(singular_values_ST_DDsU);
explained_variance_ST_DDsU_true = cumsum(singular_values_ST_DDsU_true) / sum(singular_values_ST_DDsU_true);
explained_variance_ST_DDsU_rand = cumsum(singular_values_ST_DDsU_rand) / sum(singular_values_ST_DDsU_rand);


% 累積寄与率をプロット
figure(4);
subplot(1, 2, 1)
plot(explained_variance_ST_U, 'o-', 'DisplayName', 'S3TTV');
hold on;
plot(explained_variance_ST_U_true, 'x-', 'DisplayName', 'Ground-truth');
plot(explained_variance_ST_U_rand, '^-', 'DisplayName', 'Rand');
xlabel('Number of Singular Values');
ylabel('Cumulative Explained Variance');
title('ST U');
legend;
grid on;
hold off;


subplot(1, 2, 2)
plot(explained_variance_ST_DDsU, 'o-', 'DisplayName', 'S3TTV');
hold on;
plot(explained_variance_ST_DDsU_true, 'x-', 'DisplayName', 'Ground-truth');
plot(explained_variance_ST_DDsU_rand, '^-', 'DisplayName', 'Rand');
xlabel('Number of Singular Values');
ylabel('Cumulative Explained Variance');
title('ST DDsU');
legend;
grid on;
hold off;

sgtitle('explained variance');

end