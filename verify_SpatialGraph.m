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



%% --- Guide Image Verification: Mean vs Median vs SVD ---
fprintf('\n=== Guide Image Generation Comparison ===\n');

[n1, n2, n3] = size(V);

% 1. Ground Truth Guide (Clean Mean)
% 比較基準として [0, 1] に正規化しておきます
G_clean = mean(U, 3);
G_clean = (G_clean - min(G_clean(:))) / (max(G_clean(:)) - min(G_clean(:)));

% 2. Method A: Simple Mean (Noisy)
G_mean = mean(V, 3);
G_mean = (G_mean - min(G_mean(:))) / (max(G_mean(:)) - min(G_mean(:)));

% 3. Method B: Mean + Spatial Median Filter
% ノイズ除去のために 3x3 のメディアンフィルタを適用
G_med = medfilt2(G_mean, [3 3]); 
% 正規化 (フィルタで範囲が変わる可能性があるため)
G_med = (G_med - min(G_med(:))) / (max(G_med(:)) - min(G_med(:)));

% 4. Method C: SVD (PCA) Based
% データを [画素数 x バンド数] に展開
V_flat = reshape(V, n1*n2, n3);

% SVD (Economy size) または svds (k=1) を実行
% 第1特異ベクトル（主要な構造成分）を抽出
% ※ svds は計算が速いです
[U_svd, ~, ~] = svds(V_flat, 1); 

% 画像形状に戻す
G_svd = reshape(U_svd, n1, n2);

% 【重要】SVDの符号不定性の補正
% 特異ベクトルの向きは任意なので、反転している（空が黒、地面が白など）場合があります。
% 単純な平均画像との相関を見て、負の相関なら反転させます。
if corr(G_svd(:), G_mean(:)) < 0
    G_svd = -G_svd;
end

% [0, 1] に正規化 (SVDのスコアはスケールが異なるため必須)
G_svd = (G_svd - min(G_svd(:))) / (max(G_svd(:)) - min(G_svd(:)));


%% --- Evaluation & Visualization ---

% Calculate PSNR & SSIM
psnr_mean = psnr(G_mean, G_clean);
ssim_mean = ssim(G_mean, G_clean);

psnr_med  = psnr(G_med, G_clean);
ssim_med  = ssim(G_med, G_clean);

psnr_svd  = psnr(G_svd, G_clean);
ssim_svd  = ssim(G_svd, G_clean);

fprintf('Method\t\tPSNR (dB)\tSSIM\n');
fprintf('------------------------------------\n');
fprintf('Mean\t\t%.4f\t\t%.4f\n', psnr_mean, ssim_mean);
fprintf('Mean+Med\t%.4f\t\t%.4f\n', psnr_med, ssim_med);
fprintf('SVD (PCA)\t%.4f\t\t%.4f\n', psnr_svd, ssim_svd);


% Plotting
figure('Name', 'Guide Image Comparison', 'Position', [100, 100, 1600, 800]);

% Row 1: Images
subplot(2, 4, 1); imshow(G_clean, []); title('GT Guide (Mean of Clean)');
subplot(2, 4, 2); imshow(G_mean, []);  title(sprintf('Noisy Mean\nPSNR: %.2f', psnr_mean));
subplot(2, 4, 3); imshow(G_med, []);   title(sprintf('Noisy Mean + Med\nPSNR: %.2f', psnr_med));
subplot(2, 4, 4); imshow(G_svd, []);   title(sprintf('Noisy SVD (PC1)\nPSNR: %.2f', psnr_svd));

% Row 2: Difference from GT (Error Maps)
% 差分を見やすくするために強調表示
diff_scale = 5; 

subplot(2, 4, 6); 
imagesc(abs(G_clean - G_mean) * diff_scale); axis image off; colorbar;
title('Diff: Mean'); caxis([0 1]);

subplot(2, 4, 7); 
imagesc(abs(G_clean - G_med) * diff_scale); axis image off; colorbar;
title('Diff: Mean + Med'); caxis([0 1]);

subplot(2, 4, 8); 
imagesc(abs(G_clean - G_svd) * diff_scale); axis image off; colorbar;
title('Diff: SVD'); caxis([0 1]);

colormap(gca, 'parula'); % カラーマップの設定

%% --- Optional: Zoom into a patch ---
% 細かい構造の確認用（中央付近を切り出し）
patch_size = 50;
cx = round(n2/2); cy = round(n1/2);
rect = [cx-patch_size, cy-patch_size, 2*patch_size, 2*patch_size]; % [xmin ymin w h]

G_clean_crop = imcrop(G_clean, rect);
G_mean_crop  = imcrop(G_mean, rect);
G_med_crop   = imcrop(G_med, rect);
G_svd_crop   = imcrop(G_svd, rect);

figure('Name', 'Guide Image Zoom', 'Position', [150, 150, 1200, 400]);
subplot(1,4,1); imshow(G_clean_crop, []); title('GT (Zoom)');
subplot(1,4,2); imshow(G_mean_crop, []); title('Mean (Zoom)');
subplot(1,4,3); imshow(G_med_crop, []); title('Mean+Med (Zoom)');
subplot(1,4,4); imshow(G_svd_crop, []); title('SVD (Zoom)');




%% --- Spatial Weight Accuracy: Mean vs SVD ---
% U: Clean HSI, V: Noisy HSI がロードされている前提

fprintf('\n=== Spatial Weight Quality Comparison ===\n');

% 1. Cleanな重み (理想解)
% Clean画像の平均をガイドとして作成（またはClean画像そのものから作成でもOKですが、
% Create_SpatialGraphWeight内部で平均化しているのでそのまま渡します）
[W_clean, ~] = Create_SpatialGraphWeight(U, "90");

% 2. Meanガイドから作った重み
% Create_SpatialGraphWeight は入力を平均化するので、Vをそのまま渡せば "Mean" 相当になります
[W_mean, ~] = Create_SpatialGraphWeight(V, "90");

% 3. SVDガイドから作った重み
% SVDガイド画像を作成し、それを [N1 x N2 x 1] として関数に渡します
[n1, n2, n3] = size(V);
V_flat = reshape(V, n1*n2, n3);
[U_svd, ~, ~] = svds(double(V_flat), 1);
G_svd = reshape(U_svd, n1, n2);

% SVDの符号補正 & 正規化
G_mean_ref = mean(V, 3);
if corr(G_svd(:), G_mean_ref(:)) < 0, G_svd = -G_svd; end
G_svd = (G_svd - min(G_svd(:))) / (max(G_svd(:)) - min(G_svd(:)));

% 関数に渡すため3次元に拡張 (内部で mean を取るため、同じものを複製しても結果は変わらない)
% ダミーで3バンドにして渡す
G_svd_3d = repmat(G_svd, [1 1 3]); 
[W_svd, ~] = Create_SpatialGraphWeight(G_svd_3d, "90");


% --- 評価: 重み行列の誤差比較 ---
% W は [n1, n2, n3, 4] ですが、空間グラフは全バンド共通（複製）なので
% 最初の1バンド分の重み [n1, n2, 1, 4] だけ比較すれば十分です。

W_c_slice = squeeze(W_clean(:,:,1,:)); % [n1 x n2 x 4]
W_m_slice = squeeze(W_mean(:,:,1,:));
W_s_slice = squeeze(W_svd(:,:,1,:));

% 二乗誤差 (MSE)
mse_mean = mean((W_c_slice(:) - W_m_slice(:)).^2);
mse_svd  = mean((W_c_slice(:) - W_s_slice(:)).^2);

% 相関係数 (1に近いほど構造が一致)
corr_mean = corr(W_c_slice(:), W_c_slice(:)); % Check
corr_mean_val = corr(W_m_slice(:), W_c_slice(:));
corr_svd_val  = corr(W_s_slice(:), W_c_slice(:));

fprintf('Method\t\tWeight MSE (Lower is better)\tWeight Corr (Higher is better)\n');
fprintf('------------------------------------------------------------\n');
fprintf('Mean\t\t%.6f\t\t\t%.4f\n', mse_mean, corr_mean_val);
fprintf('SVD\t\t%.6f\t\t\t%.4f\n', mse_svd, corr_svd_val);


% --- 可視化: エッジ検出能力の差 ---
% 垂直方向のエッジ重み (W(:,:,:,1)) を可視化比較
dir_idx = 1; 

figure('Name', 'Spatial Weight Comparison', 'Position', [100, 100, 1200, 600]);

subplot(1, 3, 1); 
imagesc(W_clean(:,:,1,dir_idx), [0 1]); axis image off; title('Ideal Weights (Clean)');
colorbar;

subplot(1, 3, 2); 
imagesc(W_mean(:,:,1,dir_idx), [0 1]); axis image off; 
title(sprintf('Mean Weights\nMSE: %.4f', mse_mean));
colorbar;

subplot(1, 3, 3); 
imagesc(W_svd(:,:,1,dir_idx), [0 1]); axis image off; 
title(sprintf('SVD Weights\nMSE: %.4f', mse_svd));
colorbar;

% 差分画像の表示
figure('Name', 'Weight Error Maps', 'Position', [150, 150, 1000, 500]);
subplot(1, 2, 1);
imagesc(abs(W_clean(:,:,1,dir_idx) - W_mean(:,:,1,dir_idx)), [0 1]); axis image off; colorbar;
title('Error: Mean vs Clean');

subplot(1, 2, 2);
imagesc(abs(W_clean(:,:,1,dir_idx) - W_svd(:,:,1,dir_idx)), [0 1]); axis image off; colorbar;
title('Error: SVD vs Clean');