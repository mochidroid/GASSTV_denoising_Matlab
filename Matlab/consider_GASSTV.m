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


%% Setting paramters
% sigma_sp = 0.01;
sigma_sp = "90";
sigma_l = 0.1;
num_segments = 5;


%% Creating graphs
guide_image_U = mean(U, 3);
guide_image_V = mean(V, 3);
guide_image_V_med = medfilt2(guide_image_V);


[W_spatial_U, grad_mat_2d_U, W_mat_2d_U, W_spatial_for_hist_U] = create_spatial_graph(U, sigma_sp);
[W_spatial_V, grad_mat_2d_V, W_mat_2d_V, W_spatial_for_hist_V] = create_spatial_graph(V, sigma_sp);


[W_spectral_U, segmented_image_U, segmented_spectra_U, spectra_graphs_U, spectra_diff_graphs_U] ...
    = create_spectral_graph(U, sigma_l, num_segments);
[W_spectral_V, segmented_image_V, segmented_spectra_V, spectra_graphs_V, spectra_diff_graphs_V] ...
    = create_spectral_graph(V, sigma_l, num_segments);
[W_spectral_M, segmented_image_M, segmented_spectra_M, spectra_graphs_M, spectra_diff_graphs_M] ...
    = create_spectral_graph_med(V, sigma_l, num_segments, 3);


%% Setting Operator
Dsp     = @(z) D4_Neumann_GPU(z);
Dl      = @(z) z(:, :, [2:end, end], :) - z;


%% Plotting results for spatial graph
% Showing guide images
figure('Position', [82,400,1748,900]);
subplot(3,4,1)
imshow(guide_image_U)
% imagesc(guide_image_U)
title("GT")

subplot(3,4,2)
imshow(guide_image_V)
% imagesc(guide_image_V)
title("Noisy")

subplot(3,4,3)
imagesc(abs(guide_image_U - guide_image_V))
% imagesc(guide_image_V)
title("Residual")
colorbar


% % Showing grad_mat
% grad_mat_cat_U = [grad_mat_2d_U(:,:,1) grad_mat_2d_U(:,:,2) grad_mat_2d_U(:,:,3) grad_mat_2d_U(:,:,4)];
% grad_mat_cat_V = [grad_mat_2d_V(:,:,1) grad_mat_2d_V(:,:,2) grad_mat_2d_V(:,:,3) grad_mat_2d_V(:,:,4)];
% grad_mat_cat_diff = abs(grad_mat_cat_U - grad_mat_cat_V) * 7;
% 
% 
% figure;
% subplot(2, 1, 1)
% imagesc(cat(1, grad_mat_cat_U, grad_mat_cat_V))
% colorbar
% 
% subplot(2, 1, 2)
% imagesc(grad_mat_cat_diff)
% colorbar


% Showing W_spatial
W_spatial_cat_U = [W_mat_2d_U(:,:,1) W_mat_2d_U(:,:,2) W_mat_2d_U(:,:,3) W_mat_2d_U(:,:,4)];
W_spatial_cat_V = [W_mat_2d_V(:,:,1) W_mat_2d_V(:,:,2) W_mat_2d_V(:,:,3) W_mat_2d_V(:,:,4)];


% figure('Position', [800 700 1000 500]);
subplot(3,4,[5:7, 9:11])
imagesc(cat(1, W_spatial_cat_U, W_spatial_cat_U))
title(sprintf("W_{sp} (\\sigma_{sp} = %g)", sigma_sp));
colorbar

% figure('Position', [100 100 700 500]);
subplot(3,4,[8, 12])
hold on
histogram(W_spatial_for_hist_U)
histogram(W_spatial_for_hist_V)
hold off
title(sprintf("W_{sp} (\\sigma_{sp} = %g)", sigma_sp));
legend("GT", "Noisy")


%% Plotting results for spectral graph
HSI_spectrum = reshape(U, [hsi.n1*hsi.n2, hsi.n3]);

k = 10;  % 代表スペクトルの数
[idx, Spectra_rep] = kmeans(HSI_spectrum, k, 'MaxIter', 500, 'Replicates', 5);
% C: [k, bands] 各クラスタの中心（代表スペクトル）

% 
figure('Position', [1700,500,1655,800])
subplot(2,2,1)
hold on
plot(Spectra_rep', 'k');
plot(spectra_graphs_U, 'b', 'LineWidth', 2);
plot(spectra_graphs_V, 'g', 'LineWidth', 2);
plot(spectra_graphs_M, 'r', 'LineWidth', 2);
hold off

xlabel('Band Index');
ylabel('Reflectance');
title('Representative Spectra (k-means)');
grid on;



spectra_diff_graphs_cat = cat(1, spectra_diff_graphs_U', spectra_diff_graphs_V', spectra_diff_graphs_M');


% figure;
subplot(2,2,3)
imagesc(spectra_diff_graphs_cat)
colorbar


% figure;
subplot(2,2,[2,4])
hold on
histogram(spectra_diff_graphs_U(:), 'EdgeColor', 'b', 'FaceColor', 'b')
histogram(spectra_diff_graphs_V(:), 'EdgeColor', 'g', 'FaceColor', 'g')
histogram(spectra_diff_graphs_M(:), 'EdgeColor', 'r', 'FaceColor', 'r')
hold off
title(sprintf("W_{s} (\\sigma_{s} = %g)", sigma_l));
legend("GT", "Noisy", "Med. Filt.", "Location", "west")


%% verifying med. filt. and moving average 
spectrum = spectra_graphs_V(:, 1);  % 1本の代表スペクトル

% 移動平均
s_movmean = movmean(spectrum, 5);

% メディアン
s_median = medfilt1(spectrum, 5);

% プロット比較
figure;
plot(spectrum, 'Color', [0.7 0.7 0.7]); hold on;
plot(s_movmean, 'b', 'LineWidth', 1.5);
plot(s_median, 'r', 'LineWidth', 1.5);
legend('Original', 'Moving Average', 'Median Filter');
title('Comparison of Filters');
xlabel('Band Index'); ylabel('Reflectance');
grid on;


%% Comparing Oracle and simu radius
lambda_sp_ora = sum(abs(W_spatial_U.*Dsp(U)), "all");
lambda_l_ora  = sum(abs(W_spectral_U.*Dl(U)), "all");

lambda_sp_simu = sum(abs(W_mat_2d_V.*Dsp(guide_image_V)), "all") * hsi.n3;
lambda_l_simu  = sum(abs(W_spectral_M.*Dl(segmented_spectra_M)), "all");

lambda_sp_simu_Rev = sum(abs(W_mat_2d_V.*Dsp(guide_image_V)), "all") * sum(mean(V, [1,2]));

fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n", deg.deadline_rate);

fprintf("\nRadius_sp oracle: %.4f, simu: %.4f, Rev.: %.4f\n", lambda_sp_ora, lambda_sp_simu, lambda_sp_simu_Rev)
fprintf("Radius_l oracle: %.4f, simu: %.4f\n", lambda_l_ora, lambda_l_simu)


%%
%% --- Inputs ---
% U: clean HSI      (N1 x N2 x N3)
% V: observed HSI   (N1 x N2 x N3)
% segmented_spectra_U: segment代表スペクトルを各画素に割り当てたHSI (N1 x N2 x N3)

% 例: すでにワークスペースに U, V, segmented_spectra_U がある前提

%% --- Helpers ---
% 行列カーネル（ガウシアン）から重み作成（対称化 & 対角0）
make_gaussian_weights = @(D2, epsilon) ...
    (exp(-D2 ./ (2*epsilon^2)) .* (1 - eye(size(D2))));  % D2は距離^2

% kNNでスパース化（相互kNN対称化: max で対称化）
function Wk = knn_sparsify(W, k)
    K = size(W,1);
    Wk = zeros(K);
    for i = 1:K
        [~, idx] = maxk(W(i,:), k);     % 自分以外の上位k
        Wk(i, idx) = W(i, idx);
    end
    % 対称化（maxで結合）
    Wk = max(Wk, Wk.');
    % 対角0
    Wk(1:K+1:end) = 0;
end

%% --- 1) スペクトル差分の作成（各列がノード=隣接バンド差分画像のベクトル） ---
% ノード数: P-1 (P = N3)
dU  = Dl(U);                    % N1 x N2 x N3
dV  = Dl(V);                    % N1 x N2 x N3
dS  = Dl(segmented_spectra_U);  % N1 x N2 x N3  ※グラフ構築用の差分(代表スペクトルベース)

% 最後の複製スライスはノードに使わないので 1:(P-1) を採用
P  = size(U,3);
useBands = 1:(P-1);

% 2Dに展開（各列がノード=差分ベクトル, 次元MN）
MN = size(U,1) * size(U,2);
F_clean = reshape(dU(:,:,useBands), [hsi.n1*hsi.n2, numel(useBands)]);   % from U
F_obs   = reshape(dV(:,:,useBands), [hsi.n1*hsi.n2, numel(useBands)]);   % from V
F_seg   = reshape(dS(:,:,useBands), [hsi.n1*hsi.n2, numel(useBands)]);   % for graph weights (segmented_spectra_U)

%% --- 2) 重み行列 W の計算（代表スペクトルベース差分から） ---
% ノード間ユークリッド距離の二乗行列 D2 を計算
% ここでは列間距離: D2(i,j) = ||F_seg(:,i) - F_seg(:,j)||_2^2
G = F_seg;                        % (MN x K), K = P-1
K = size(G,2);
GTG = G.'*G;                      % (K x K)
nrm2 = diag(GTG);                 % 各列の二乗ノルム
D2 = (nrm2 + nrm2.') - 2*GTG;     % 距離^2の行列（対称, 対角0）

% epsilonの自動設定（ゼロ以外の中央値を採用）
d2_vec = D2(triu(true(K),1));  % 上三角のオフ対角要素
d2_vec = d2_vec(d2_vec>0);
if isempty(d2_vec), d2_vec = 1; end
epsilon = sqrt(median(d2_vec));   % epsilon ~ 距離の中央値の平方根

% ガウシアンカーネルで重み作成
W_full = make_gaussian_weights(D2, epsilon);  % 密行列

% スパース化（相互kNN）
k = min(10, K-1);                 % デフォルト: 上位10
W = knn_sparsify(W_full, k);      % 最終的に使う重み（スパース）

% グラフラプラシアン
d = sum(W,2);
L = diag(d) - W;

%% --- 3) GLR評価関数の計算 ---
% SG(F) = Tr(F L F^T) = sum_i F(:,i)^T * (sum_j W_ij [F(:,i)-F(:,j)])
% ここでは clean と observed で比較
SG_clean = trace(F_clean * L * F_clean.');
SG_obs   = trace(F_obs   * L * F_obs.');

% 正規化（ノード数や重み総和の影響を緩和したい場合）
Wsum = sum(W(:))/2;     % 無向グラフの総重み
if Wsum == 0, Wsum = 1; end
SG_clean_norm = SG_clean / Wsum;
SG_obs_norm   = SG_obs   / Wsum;

%% --- 4) 参考: 指標の表示 ---
fprintf('Nodes (K=P-1): %d\n', K);
fprintf('Nonzeros in W: %d (density %.3f)\n', nnz(W), nnz(W)/numel(W));
fprintf('epsilon (auto): %.4g\n', epsilon);
fprintf('SG_clean: %.6g  | SG_obs: %.6g\n', SG_clean, SG_obs);
fprintf('SG_clean_norm: %.6g  | SG_obs_norm: %.6g\n', SG_clean_norm, SG_obs_norm);


%% ====== 3) 可視化（ヒートマップとグラフ図） ======
% (a) 重み行列のヒートマップ
figure; imagesc(W); axis image; colorbar;
xlabel('ΔBand index'); ylabel('ΔBand index');
title('Spectral Difference Graph (G1): weight matrix W');
set(gca, 'XTick', 1:K, 'YTick', 1:K);

% (b) バンド差分インデックスを付けたグラフ図（円形）
G = graph(W, 'upper');
nodeNames = arrayfun(@(i) sprintf('\\Delta(%d,%d)', i, i+1), 1:K, 'UniformOutput', false);
figure;
h = plot(G, 'Layout', 'circle', 'NodeLabel', nodeNames);
maxW = max(W(:));
if maxW > 0
    h.LineWidth = 5 * (G.Edges.Weight / maxW);
    h.EdgeCData = G.Edges.Weight; colormap(jet); colorbar;
end
title('G1 graph (edge width/color = weight)');


%% Creating spatial graph
function [W_spatial, grad_mat_2d, W_mat_2d, W_spatial_for_hist] = create_spatial_graph(X, sigma_sp)
guide_image = single(mean(X, 3));
[n1, n2, n3] = size(X);


diff_1_0 = cat(1, guide_image(2:n1, :), Inf(1, n2));
diff_0_1 = cat(2, guide_image(:, 2:n2), Inf(n1, 1));

diff_1_1 = cat(2, diff_1_0(:, 2:n2), Inf(n1, 1));
diff_1_2 = cat(1, Inf(1, n2), diff_0_1(1:n1-1, :));

diff_mat = [diff_1_0(:), diff_0_1(:), diff_1_1(:), diff_1_2(:)];
grad_mat = guide_image(:).*ones(1, 4) - diff_mat;

grad_mat_2d = reshape(grad_mat, [n1, n2, 4]);                         

% dist_mat = [1, 1, 2, 2];
if sigma_sp == "med"
    diff_vals = abs(grad_mat(~isinf(grad_mat)));
    val_sigma_sp = median(diff_vals(:));
elseif sigma_sp == "90"
    finite_mask = isfinite(grad_mat);
    g = abs(grad_mat(finite_mask));
    g0 = prctile(g, 90);   % representative gradient magnitude
    w0 = 0.2;              % desired weight at g0
    val_sigma_sp = max( single( g0 / sqrt(2*log(1/w0)) ), eps('single') );
else
    val_sigma_sp = sigma_sp;
end


% W_mat_tmp = exp(-(grad_mat.^2)./(sigma_x.^2)/2).*exp(-dist_mat./(sigma_l.^2)/2);
W_mat_tmp = exp(-(grad_mat.^2)/(val_sigma_sp^2)/2);
% W_mat_tmp = exp(-abs(grad_mat)/(sigma_x));
W_mat = repmat(W_mat_tmp, n3, 1);

W_mat_2d = reshape(W_mat_tmp, [n1, n2, 4]);

W_spatial_for_hist = {W_mat_2d(1:end-1, :, 1), ...
                        W_mat_2d(:, 1:end-1, 2), ...
                        W_mat_2d(1:end-1, 1:end-1, 3), ...
                        W_mat_2d(2:end, 1:end-1, 4)};
W_spatial_for_hist = [W_spatial_for_hist{1}(:); W_spatial_for_hist{2}(:); W_spatial_for_hist{3}(:); W_spatial_for_hist{4}(:)];

W_spatial = reshape(W_mat, [n1, n2, n3, 4]);
W_spatial = single(W_spatial);


fprintf('sigma_sp = %g\n', val_sigma_sp)

end


%% Creating spectral graph
function [W_spectral, segmented_image, segmented_spectra, spectra_graphs, spectra_diff_graphs] ...
    = create_spectral_graph(X, sigma_s, num_segments)
guide_image = single(mean(X, 3));
[n1, n2, n3] = size(X);

segmented_image = imsegkmeans(guide_image, num_segments);
segmented_image = reshape(segmented_image, size(guide_image));

W_spectral = zeros([n1, n2, n3]);
segmented_spectra = zeros([n1, n2, n3]);
spectra_graphs = zeros([n3, num_segments]);
spectra_diff_graphs = zeros([n3, num_segments]);

for i = 1:num_segments
    % セグメントごとのピクセルインデックスを抽出
    mask = (segmented_image == i);
    mask3D = repmat(mask, [1,1,n3]);
    % セグメントに対応するHSデータを空間方向に平均
    HSI_segmented = X.*mask3D;
    % spectral_graph = mean(HSI_segmented, [1,2]);  % 各バンドの平均
    spectral_graph = sum(HSI_segmented, [1,2]) / sum(mask, "all");  % 各バンドの平均
    spectra_graphs(:, i) = spectral_graph;

    spectral_diff = spectral_graph(:, :,[2:end, end]) - spectral_graph;
    spectral_diff_graph = exp(-(spectral_diff.^2)/(sigma_s^2)/2);
    spectra_diff_graphs(:, i) = spectral_diff_graph;


    W_spectral = W_spectral + mask3D.*spectral_diff_graph;
    segmented_spectra = segmented_spectra + mask3D.*spectral_graph;
end

W_spectral = single(W_spectral);
segmented_spectra = single(segmented_spectra);

end


%% Creating spectral graph med. filt.
function [W_spectral, segmented_image, segmented_spectra, spectra_graphs, spectra_diff_graphs] ...
    = create_spectral_graph_med(X, sigma_s, num_segments, med_order)
guide_image = single(mean(X, 3));
[n1, n2, n3] = size(X);

segmented_image = imsegkmeans(guide_image, num_segments);
segmented_image = reshape(segmented_image, size(guide_image));

W_spectral = zeros([n1, n2, n3]);
segmented_spectra = zeros([n1, n2, n3]);
spectra_graphs = zeros([n3, num_segments]);
spectra_diff_graphs = zeros([n3, num_segments]);

for i = 1:num_segments
    % セグメントごとのピクセルインデックスを抽出
    mask = (segmented_image == i);
    mask3D = repmat(mask, [1,1,n3]);
    % セグメントに対応するHSデータを空間方向に平均
    HSI_segmented = X.*mask3D;
    % spectral_graph = mean(HSI_segmented, [1,2]);  % 各バンドの平均
    spectral_graph = sum(HSI_segmented, [1,2]) / sum(mask, "all");  % 各バンドの平均
    spectral_graph = medfilt1(spectral_graph, med_order);
    spectra_graphs(:, i) = spectral_graph;

    spectral_diff = spectral_graph(:, :,[2:end, end]) - spectral_graph;
    spectral_diff_graph = exp(-(spectral_diff.^2)/(sigma_s^2)/2);
    spectra_diff_graphs(:, i) = spectral_diff_graph;


    W_spectral = W_spectral + mask3D.*spectral_diff_graph;
    segmented_spectra = segmented_spectra + mask3D.*spectral_graph;
end

W_spectral = single(W_spectral);
segmented_spectra = single(segmented_spectra);

end