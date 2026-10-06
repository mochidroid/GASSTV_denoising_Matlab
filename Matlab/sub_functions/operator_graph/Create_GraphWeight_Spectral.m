function [W_spectral] = Create_GraphWeight_Spectral(X, sigma_sp, num_segments)
[n1, n2, n3] = size(X);

guide_image = mean(X, 3);
% segmented_image = kmeans(guide_image(:), num_segments);
segmented_image = imsegkmeans(guide_image, num_segments);
segmented_image = reshape(segmented_image, size(guide_image));

W_spectral = zeros([n1, n2, n3]);

for i = 1:num_segments
    % セグメントごとのピクセルインデックスを抽出
    mask = (segmented_image == i);
    mask3D = repmat(mask, [1,1,n3]);
    % セグメントに対応するHSデータを空間方向に平均
    HSI_segmented = X.*mask3D;
    % spectral_graph = mean(HSI_segmented, [1,2]);  % 各バンドの平均
    spectral_graph = sum(HSI_segmented, [1,2]) / sum(mask, "all");  % 各バンドの平均
    spectral_diff = spectral_graph(:, :,[2:end, end]) - spectral_graph;
    spectral_diff_graph = exp(-(spectral_diff.^2)/(2*sigma_sp^2));
    W_spectral = W_spectral + mask3D.*spectral_diff_graph;
end

W_spectral = single(W_spectral);