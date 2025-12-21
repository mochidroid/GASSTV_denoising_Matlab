close all
addpath(genpath("./sub_functions"));


%% Selecting conditions
images = {...
    "IndianPines", ...
    "Suwannee", ...
};

idx_image = 1;
% idx_image = 2;


%% Generating observation
image = images{idx_image};

[V, hsi] = Load_real_HSI(image);
V = single(V);

fprintf('Running Verification on: %s\n', image);


%% Setting Parameters
sigma_sp_opts = "90";   % 空間グラフ用 sigma ("med" or "90" or numeric)
num_segments = single(5);
k_lap        = 30;
sigma_l      = "med"; % スペクトルグラフ用
order_filt   = single(3);     % 代表スペクトルのメディアンフィルタサイズ


%% 2. Spatial Graph Comparison
fprintf('\n--- Calculating Spatial Graphs ---\n');

% Noisy画像から作成
[Wsp_noisy, rho_noisy] = Create_SpatialGraphWeight(V, sigma_sp_opts);

% ガイド画像の取得（Create_SpatialGraphWeight内で平均をとっているため、ここでも再計算）
guide_V = mean(V, 3);


% 重みの差分 (Wspは [n1, n2, n3, 4] だが、n3方向は複製なので1スライスだけ見る)
Wsp_noisy_slice = squeeze(Wsp_noisy(:,:,1,:));

%% 3. Spectral Graph (GLR) Comparison
fprintf('--- Calculating Spectral Graphs (GLR) ---\n');

% Noisy
[L_noisy, ~, ~, info_noisy] = ...
    Create_SpectralGraphLaplacian(V, num_segments, sigma_l, k_lap, order_filt);

% 代表スペクトルの取得
R_noisy = info_noisy.representative;

% 注: ノイズによりkmeansの結果（ラベル割り当て）が変わるため、
% 単純なセグメントID順の比較はずれる可能性があります。
% ここでは、「最も画素数が多いセグメント（背景など）」同士を比較してみます。
[~, idx_max_noisy] = max(info_noisy.segment_sizes);

% そのセグメントのグラフラプラシアン L (K x K)
Ln_seg = L_noisy(:,:,idx_max_noisy);

% 重み行列 W の復元 (W = diag(L) - L, off-diagonal only)
Wn_seg = diag(diag(Ln_seg)) - Ln_seg;

%% 4. Visualization (Real Data Analysis)
% Clean画像がないため、Noisyデータの統計と構造のみを可視化します

% === Figure 1: Spatial Graph Analysis ===
figure('Name', 'Spatial Graph Analysis (Real Data)', 'Position', [100, 100, 1200, 500]);

% 1. ガイド画像 (Mean Image)
subplot(1, 3, 1); 
imshow(guide_V, []); 
title('Guide Image (Mean of V)');
colorbar;

% 2. 空間重みの可視化 (例: 水平方向の重み channel 2)
w_idx = 2; 
subplot(1, 3, 2); 
imagesc(Wsp_noisy_slice(:,:,w_idx), [0 1]); 
colorbar; 
title(sprintf('Spatial Weights W_{sp} (Dir %d)', w_idx)); 
axis image off;

% 3. 重みヒストグラム
subplot(1, 3, 3); 
histogram(Wsp_noisy_slice(:), 100, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
grid on; 
title('Histogram of Spatial Weights');
xlabel('Weight Value'); ylabel('PDF');
xlim([0 1]);


% === Figure 2: Spectral Graph (GLR) Analysis ===
figure('Name', 'Spectral Graph Analysis (Real Data)', 'Position', [150, 150, 1000, 500]);
col_n = 'r';  % Noisy color

% 1. 代表スペクトルのプロット
subplot(1, 2, 1);
% 全セグメントのスペクトルを薄く表示
plot(R_noisy', '-', 'Color', [1 0 0 0.2]); hold on;
% 最大セグメントを太く強調
plot(R_noisy(idx_max_noisy, :), 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Max Segment');
% 平均スペクトル
plot(mean(R_noisy,1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Global Mean');
legend('Location', 'best');
grid on; axis tight;
title('Representative Spectra (GLR)');
xlabel('Band Index'); ylabel('Reflectance');

% 2. 最大セグメントの重み行列 (Adjacency Matrix)
subplot(1, 2, 2);
imagesc(Wn_seg); 
colorbar; 
axis square; 
title(sprintf('Adjacency W (Seg %d: Largest)', idx_max_noisy)); 
xlabel('Band Idx'); ylabel('Band Idx');


%% 5. Verification Print Output
fprintf('\n--- Verification Results (Real Data) ---\n');
% Cleanとの比較系は削除し、単体の統計量を表示
fprintf('Spatial Radius (rho): %.4e\n', rho_noisy);
fprintf('Max Segment Index: %d (Size: %d pixels)\n', idx_max_noisy, info_noisy.segment_sizes(idx_max_noisy));
fprintf('Max Eigenvalue (L) of Max Seg: %.4f\n', max(real(eig(Ln_seg))));


%% --- 6. Visualization: Per-Segment Stacked Analysis ---
% 表示するセグメントの範囲
disp_segments = 1:min(10, size(L_noisy, 3)); 
num_disp = length(disp_segments);

figure('Name', 'Per-Segment Analysis', 'Position', [50, 50, 800, 150*num_disp]);
t = tiledlayout(num_disp, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
c_lim = [0 1]; % 重みの表示レンジ

for i = 1:num_disp
    s = disp_segments(i);
    
    % --- Data Extraction ---
    % 1. Spectrum
    spec_n = info_noisy.representative(s, :);
    
    % 2. Weights (Reconstruct W from L)
    Ln = L_noisy(:,:,s);
    Wn = diag(diag(Ln)) - Ln;
    Wn(1:size(Wn,1)+1:end) = 0; % Remove numerical noise on diagonal
    
    % --- Plotting ---
    
    % Col 1: Representative Spectrum
    nexttile(t);
    plot(spec_n, 'r-', 'LineWidth', 1.5);
    grid on; axis tight;
    title(sprintf('Seg %d Spectrum', s), 'FontSize', 9);
    % Y軸ラベルの代わりにセグメント番号を表示
    ylabel(sprintf('Seg %d', s), 'FontWeight', 'bold');
    
    % Col 2: Weight Matrix
    nexttile(t);
    imagesc(Wn, c_lim); 
    colorbar;
    axis square; axis off;
    title(sprintf('Seg %d Weights', s), 'FontSize', 9);
    
    % 統計情報の表示
    text(0, size(Wn,1)+5, ...
        sprintf('Mean W: %.2f', mean(Wn(:))), ...
        'FontSize', 8, 'Interpreter', 'none');
end

xlabel(t, 'Band Index');