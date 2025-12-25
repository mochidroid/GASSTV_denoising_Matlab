function [Wsp, rho_sp] = Create_SpatialGraphWeight(X, sigma_sp)
fprintf("~ Creating spatial graph ~\n")
% Creating graph based weight matrix
guide_image = mean(X, 3);
[n1,n2,n3] = size(X);

diff_1_0 = cat(1, guide_image(2:n1, :), Inf(1, n2));
diff_0_1 = cat(2, guide_image(:, 2:n2), Inf(n1, 1));

diff_1_1 = cat(2, diff_1_0(:, 2:n2), Inf(n1, 1));
diff_1_2 = cat(1, Inf(1, n2), diff_0_1(1:n1-1, :));

diff_mat = [diff_1_0(:), diff_0_1(:), diff_1_1(:), diff_1_2(:)];
grad_mat = guide_image(:).*ones(1, 4) - diff_mat;

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

fprintf("sigma sp: %5f\n", val_sigma_sp);

% W_mat_tmp = exp(-(grad_mat.^2)./(sigma_x.^2)/2).*exp(-dist_mat./(sigma_l.^2)/2);
W_mat_tmp = exp(-(grad_mat.^2)/(val_sigma_sp^2)/2);
% W_mat_tmp = exp(-abs(grad_mat)/(sigma_x));
W_mat = repmat(W_mat_tmp, n3, 1);

Wsp = reshape(W_mat, [n1, n2, n3, 4]);
Wsp = single(Wsp);

% implay(cat(2, W(:,:,:,1), W(:,:,:,2), W(:,:,:,3), W(:,:,:,4)))

% Estimating radius of spatial graph
% (A) ガイド画像の空間グラフ差分値 (1バンドあたり)
%     Inf（境界）を除外して、重み付き差分の総和を計算
valid_mask = isfinite(grad_mat);
guide_weighted_diff = W_mat_tmp(valid_mask) .* abs(grad_mat(valid_mask));
val_sp_guide = sum(guide_weighted_diff(:));

% (B) 各バンド毎の平均輝度値ベクトルとその総和
%     spatial mean per band -> sum over bands
mean_lum_vec = squeeze(mean(mean(X, 1), 2)); 
sum_lum = sum(mean_lum_vec);

rho_sp = single(val_sp_guide * sum_lum);

end