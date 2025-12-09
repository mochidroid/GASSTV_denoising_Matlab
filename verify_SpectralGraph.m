close all
addpath(genpath("./sub_functions"));


%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % g0.05 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
    {0.1,   0,     0,     0,    0.01},  ... % g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % g0.1 ps0.05 pt0.05 pd0.01
};

idx_noise_condition = 7;

images = {... 
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

idx_image = 1;


%% Generating observation
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};

[U, hsi] = Load_HSI(image);
noise_seed = "default";
[V, deg] = Generate_obsv_for_denoising(U, deg, noise_seed);

U = single(U);
V = single(V);

fprintf('Running Verification on: %s\n', image);
fprintf('Noise: Gaussian=%.2f, Sparse=%.2f, Stripe=%.2f\n', ...
    deg.gaussian_sigma, deg.sparse_rate, deg.stripe_rate);


%% --- Segmentation Comparison: K-means vs Superpixels (SLIC) ---
% 既に U (Clean), V (Noisy) がワークスペースにある前提
% もしまだならロードしてください
% [U, hsi] = Load_HSI("JasperRidge"); 
% [V, deg] = Generate_obsv_for_denoising(U, deg, "default");

fprintf('\n=== Segmentation Method Comparison ===\n');

% パラメータ設定
target_num_segments = 5; % 比較のため少し多めに設定
med_filt_order = 5;        % 以前の検証で決まった平滑化強度

% ガイド画像 (Mean) の作成
Guide_V = mean(V, 3);
% SLIC用に [0,1] 正規化 (必須ではないが推奨)
Guide_V = (Guide_V - min(Guide_V(:))) / (max(Guide_V(:)) - min(Guide_V(:)));

%% 1. Method A: K-means (imsegkmeans)
tic;
L_km = imsegkmeans(Guide_V, target_num_segments);
time_km = toc;
num_km = max(L_km(:));

%% 2. Method B: Superpixels (SLIC)
% MATLABのsuperpixels関数を使用
tic;
[L_sp, num_sp] = superpixels(Guide_V, target_num_segments);
time_sp = toc;

fprintf('Method\t\tSegments\tTime(s)\n');
fprintf('--------------------------------\n');
fprintf('K-means\t\t%d\t\t%.4f\n', num_km, time_km);
fprintf('SLIC\t\t%d\t\t%.4f\n', num_sp, time_sp);


%% 3. Evaluate Representative Spectra Accuracy
% 各セグメントについて、「Noisy画像からの中央値」と「Clean画像からの中央値」を比較

% --- K-means Evaluation ---
[err_km, homo_km] = eval_segmentation(U, V, L_km, num_km, med_filt_order);

% --- SLIC Evaluation ---
[err_sp, homo_sp] = eval_segmentation(U, V, L_sp, num_sp, med_filt_order);

fprintf('\nEvaluation Results (Lower is better):\n');
fprintf('Method\t\tSpec. MSE\tHomogeneity (Std)\n');
fprintf('--------------------------------------------\n');
fprintf('K-means\t\t%.6f\t%.6f\n', mean(err_km), mean(homo_km));
fprintf('SLIC\t\t%.6f\t%.6f\n', mean(err_sp), mean(homo_sp));


%% 4. Visualization

figure('Name', 'Segmentation Comparison', 'Position', [100, 100, 1200, 800]);

% Row 1: Label Maps
subplot(2, 2, 1);
imshow(label2rgb(L_km, 'jet', 'k', 'shuffle'));
title(sprintf('K-means Labels (N=%d)', num_km));

subplot(2, 2, 2);
imshow(label2rgb(L_sp, 'jet', 'k', 'shuffle'));
title(sprintf('Superpixels (SLIC) Labels (N=%d)', num_sp));

% Row 2: Error Distribution (Histogram)
subplot(2, 2, 3);
h1 = histogram(err_km, 20, 'Normalization', 'pdf'); hold on;
h2 = histogram(err_sp, 20, 'Normalization', 'pdf');
h1.FaceAlpha = 0.5; h2.FaceAlpha = 0.5;
legend('K-means', 'SLIC');
title('Spectral Error Distribution (Noisy vs Clean)');
xlabel('MSE of Representative Spectrum');

% Row 2: Homogeneity (Clean image std dev inside segments)
subplot(2, 2, 4);
% 各セグメントの「中身の均一性」を表示
b = bar([mean(homo_km), mean(homo_sp)]);
xticklabels({'K-means', 'SLIC'});
title('Intra-segment Std Dev (Lower = More Homogeneous)');
ylabel('Average Std Dev (Clean)');


%% --- Helper Function ---
function [mse_list, homo_list] = eval_segmentation(Clean, Noisy, LabelMap, num_seg, filter_order)
    [~, ~, n3] = size(Clean);
    
    % 高速化のために2次元配列化
    C_flat = reshape(Clean, [], n3);
    N_flat = reshape(Noisy, [], n3);
    L_flat = LabelMap(:);
    
    mse_list = zeros(num_seg, 1);
    homo_list = zeros(num_seg, 1);
    
    for i = 1:num_seg
        idx = (L_flat == i);
        if ~any(idx), continue; end
        
        % 1. Clean画像におけるセグメント内の均一性 (標準偏差の平均)
        % 物理的に同じ物質を選べているかの指標
        pixels_c = C_flat(idx, :);
        std_c = std(pixels_c, 0, 1);
        homo_list(i) = mean(std_c);
        
        % 2. 代表スペクトルの再現誤差
        % Clean (Ground Truth) Representative
        rep_c = median(pixels_c, 1);
        
        % Noisy Representative (Proposed Method)
        pixels_n = N_flat(idx, :);
        rep_n = median(pixels_n, 1);
        % Apply smoothing (optional, based on your pipeline)
        if filter_order > 1
            rep_n = medfilt1(rep_n, filter_order);
        end
        
        % MSE Calculation
        mse_list(i) = mean((rep_c - rep_n).^2);
    end
end