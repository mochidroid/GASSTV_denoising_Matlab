function [W] = Calc_Weight(HSI_noisy, params)
boundary    = params.boundary;
sigma_x     = params.sigma_x;

[~, ~, n3]  = size(HSI_noisy);

%% Generating guide map
guide_map = mean(HSI_noisy, 3);


%% Calculating Diff
switch boundary
    case 'Neumann'
        Diff_v = guide_map([2:end, end], :) - guide_map;
        Diff_h = guide_map(:, [2:end, end]) - guide_map;

    case 'circulant'
        Diff_v = guide_map([2:end, 1], :) - guide_map;
        Diff_h = guide_map(:, [2:end, 1]) - guide_map;
end

%% Generating Weight matrix
W_v = exp(-(abs(Diff_v)/sigma_x));
W_h = exp(-(abs(Diff_h)/sigma_x));

W = cat(4, repmat(W_v, [1,1,n3]), repmat(W_h, [1,1,n3]));

%% Confirming Weight matrix
% figure;
% hist(W(:))
% figure;
% imagesc(cat(2, W_v, W_h))

