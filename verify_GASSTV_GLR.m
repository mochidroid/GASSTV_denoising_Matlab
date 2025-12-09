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


%% Setting Parameters
sigma_sp_opts = "90";   % 空間グラフ用 sigma ("med" or "90" or numeric)
opts_GLR.num_segments = single(5);
opts_GLR.k_lap        = 30;
opts_GLR.sigma_l      = "med"; % スペクトルグラフ用
opts_GLR.order_filt   = single(3);     % 代表スペクトルのメディアンフィルタサイズ


%% 2. Spatial Graph Comparison
fprintf('\n--- Calculating Spatial Graphs ---\n');

% Clean画像から作成
[Wsp_clean, rho_clean] = Create_SpatialGraphWeight(U, sigma_sp_opts);
% Noisy画像から作成
[Wsp_noisy, rho_noisy] = Create_SpatialGraphWeight(V, sigma_sp_opts);

% ガイド画像の取得（Create_SpatialGraphWeight内で平均をとっているため、ここでも再計算）
guide_U = mean(U, 3);
guide_V = mean(V, 3);

% メトリクス計算
diff_guide = abs(guide_U - guide_V);
psnr_guide = 10 * log10(1 / mean(diff_guide(:).^2));

% 重みの差分 (Wspは [n1, n2, n3, 4] だが、n3方向は複製なので1スライスだけ見る)
Wsp_clean_slice = squeeze(Wsp_clean(:,:,1,:)); % [n1, n2, 4]
Wsp_noisy_slice = squeeze(Wsp_noisy(:,:,1,:));
diff_Wsp = abs(Wsp_clean_slice - Wsp_noisy_slice);

%% 3. Spectral Graph (GLR) Comparison
fprintf('--- Calculating Spectral Graphs (GLR) ---\n');

% Clean
[L_clean, ~, ~, info_clean] = Create_SpectralGraphLaplacian(U, opts_GLR);
% Noisy
[L_noisy, ~, ~, info_noisy] = Create_SpectralGraphLaplacian(V, opts_GLR);

% 代表スペクトルの取得
R_clean = info_clean.representative; % [S x K] (or S x N3 depending on func)
R_noisy = info_noisy.representative;

% 注: ノイズによりkmeansの結果（ラベル割り当て）が変わるため、
% 単純なセグメントID順の比較はずれる可能性があります。
% ここでは、「最も画素数が多いセグメント（背景など）」同士を比較してみます。
[~, idx_max_clean] = max(info_clean.segment_sizes);
[~, idx_max_noisy] = max(info_noisy.segment_sizes);

% そのセグメントのグラフラプラシアン L (K x K)
Lc_seg = L_clean(:,:,idx_max_clean);
Ln_seg = L_noisy(:,:,idx_max_noisy);

% 重み行列 W の復元 (W = diag(L) - L, off-diagonal only)
Wc_seg = diag(diag(Lc_seg)) - Lc_seg;
Wn_seg = diag(diag(Ln_seg)) - Ln_seg;

%% 4. Visualization

% === Figure 1: Spatial Graph Analysis ===
figure('Name', 'Spatial Graph Analysis', 'Position', [100, 100, 1600, 900]);

% 1. ガイド画像比較
subplot(2, 4, 1); imshow(guide_U, []); title('Guide Image (Clean)');
subplot(2, 4, 2); imshow(guide_V, []); title(sprintf('Guide Image (Noisy)\nPSNR: %.2f dB', psnr_guide));
subplot(2, 4, 3); imagesc(diff_guide); colorbar; title('Guide Difference'); axis image off;

% 2. 空間重みの比較 (例: 水平方向の重み channel 2)
% ※ channel定義は実装依存ですが、通常 1:Vertical, 2:Horizontal と仮定
w_idx = 2; 
subplot(2, 4, 5); imagesc(Wsp_clean_slice(:,:,w_idx), [0 1]); colorbar; 
title(sprintf('W_{sp} Clean (Dir %d)', w_idx)); axis image off;

subplot(2, 4, 6); imagesc(Wsp_noisy_slice(:,:,w_idx), [0 1]); colorbar; 
title(sprintf('W_{sp} Noisy (Dir %d)', w_idx)); axis image off;

subplot(2, 4, 7); imagesc(diff_Wsp(:,:,w_idx)); colorbar; 
title('Weight Difference'); axis image off;

% 3. 重みヒストグラム
subplot(2, 4, [4, 8]); 
hold on;
histogram(Wsp_clean_slice(:), 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Clean');
histogram(Wsp_noisy_slice(:), 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Noisy');
hold off;
legend; grid on; title('Histogram of Spatial Weights');
xlabel('Weight Value'); ylabel('PDF');


% === Figure 2: Spectral Graph (GLR) Analysis ===
figure('Name', 'Spectral Graph Analysis', 'Position', [150, 150, 1400, 800]);

% プロット設定：見やすくするために線を太く、色を明るく
lw = 1.5;     % LineWidth
col_c = 'c';  % Clean color (Cyan) -> 黒背景で見やすい
col_n = 'r';  % Noisy color (Red)

% 1. 代表スペクトルの比較
subplot(2, 3, 1);
% 薄く全線をプロット
plot(R_clean', '-', 'Color', [0 1 1 0.3]); hold on; 
% 平均を太くプロット
plot(mean(R_clean,1), '-', 'Color', col_c, 'LineWidth', 2, 'DisplayName', 'Clean Mean'); 
title('Representative Spectra (Clean)'); grid on; axis tight;
set(gca, 'Color', [0.1 0.1 0.1]); % 背景を少し明るく

subplot(2, 3, 2);
plot(R_noisy', '-', 'Color', [1 0 0 0.3]); hold on;
plot(mean(R_noisy,1), '-', 'Color', col_n, 'LineWidth', 2, 'DisplayName', 'Noisy Mean');
title('Representative Spectra (Noisy)'); grid on; axis tight;
set(gca, 'Color', [0.1 0.1 0.1]);

subplot(2, 3, 3);
% 特定セグメント（最大領域）のスペクトル比較
plot(R_clean(idx_max_clean, :), 'Color', col_c, 'LineWidth', 2, 'DisplayName', 'Clean (Max Seg)');
hold on;
plot(R_noisy(idx_max_noisy, :), 'Color', col_n, 'LineWidth', 2, 'DisplayName', 'Noisy (Max Seg)');
legend('Location', 'best', 'TextColor', 'w'); 
title('Comparison: Largest Segment Spectrum'); grid on; axis tight;
set(gca, 'Color', [0.1 0.1 0.1]);


% 2. スペクトルグラフ（重み行列）の比較
% 最大セグメントの重み行列 W (K x K) を表示
subplot(2, 3, 4);
imagesc(Wc_seg); colorbar; axis square; 
title('Adjacency W (Clean, Max Seg)'); xlabel('Band Idx'); ylabel('Band Idx');

subplot(2, 3, 5);
imagesc(Wn_seg); colorbar; axis square; 
title('Adjacency W (Noisy, Max Seg)'); xlabel('Band Idx'); ylabel('Band Idx');

subplot(2, 3, 6);
% 重みの散布図比較（構造がどれくらい保たれているか）
% 行列サイズが同じ場合のみ
if size(Wc_seg) == size(Wn_seg)
    scatter(Wc_seg(:), Wn_seg(:), 10, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
    hold on; plot([0 1], [0 1], 'r--'); % y=x line
    axis([0 1 0 1]); grid on;
    xlabel('Clean Weights'); ylabel('Noisy Weights');
    title('Weight Correlation');
else
    text(0.5, 0.5, 'Size mismatch due to different K', 'HorizontalAlignment', 'center');
end

%% 5. Verification Print Output
fprintf('\n--- Verification Results ---\n');
fprintf('Guide Image PSNR: %.2f dB\n', psnr_guide);
fprintf('Spatial Radius (rho) - Clean: %.4e, Noisy: %.4e (Ratio: %.2f)\n', ...
    rho_clean, rho_noisy, rho_noisy/rho_clean);
fprintf('Max Segment Index - Clean: %d, Noisy: %d\n', idx_max_clean, idx_max_noisy);
fprintf('Max Eigenvalue (L) - Clean (Max Seg): %.4f\n', max(real(eig(Lc_seg))));
fprintf('Max Eigenvalue (L) - Noisy (Max Seg): %.4f\n', max(real(eig(Ln_seg))));



%% --- 6. Visualization: Per-Segment Stacked Analysis ---
% 表示するセグメントの範囲（全部だと重くなるので適宜調整してください）
% 例: 1:5 (最初の5個), または 1:S_seg (全部)
disp_segments = 1:min(10, size(L_clean, 3)); 
num_disp = length(disp_segments);

figure('Name', 'Per-Segment Analysis', 'Position', [50, 50, 1200, 150*num_disp]);
t = tiledlayout(num_disp, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

% 共通設定
c_lim = [0 1]; % 重みの表示レンジ [0, 1]

for i = 1:num_disp
    s = disp_segments(i);
    
    % --- Data Extraction ---
    % 1. Spectra
    spec_c = info_clean.representative(s, :);
    spec_n = info_noisy.representative(s, :);
    
    % 2. Weights (Reconstruct W from L)
    % L = D - W => W = diag(diag(L)) - L (off-diagonal)
    Lc = L_clean(:,:,s);
    Wc = diag(diag(Lc)) - Lc;
    % 数値誤差で対角成分が残ることがあるのでゼロにする
    Wc(1:size(Wc,1)+1:end) = 0; 
    
    Ln = L_noisy(:,:,s);
    Wn = diag(diag(Ln)) - Ln;
    Wn(1:size(Wn,1)+1:end) = 0;
    
    % --- Plotting ---
    
    % Col 1: Representative Spectra Comparison
    nexttile(t);
    plot(spec_c, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Clean'); hold on;
    plot(spec_n, 'r-', 'LineWidth', 1.0, 'DisplayName', 'Noisy');
    grid on; axis tight;
    % Y軸の範囲を合わせる（平坦かどうかの確認用）
    ylim([min([spec_c(:); spec_n(:)])-0.05, max([spec_c(:); spec_n(:)])+0.05]);
    if i == 1
        title('Rep. Spectra'); 
        legend('Location','best');
    end
    ylabel(sprintf('Seg %d', s), 'FontWeight', 'bold');
    
    % Col 2: Clean Weight Matrix
    nexttile(t);
    imagesc(Wc, c_lim); 
    axis square; axis off;
    if i == 1, title('Clean Weights W'); end
    % colorbar; % スペース節約のため省略（色は黄色=1, 青=0）
    
    % Col 3: Noisy Weight Matrix
    nexttile(t);
    imagesc(Wn, c_lim); 
    axis square; axis off;
    if i == 1, title('Noisy Weights W'); end
    
    % デバッグ情報: 重みの最大・最小・平均を表示
    text(0, size(Wn,1)+5, ...
        sprintf('Mean W: %.2f', mean(Wn(:))), ...
        'FontSize', 8, 'Interpreter', 'none');
end

% colormap(t, 'parula'); % 重み用のカラーマップ
xlabel(t, 'Band Index');