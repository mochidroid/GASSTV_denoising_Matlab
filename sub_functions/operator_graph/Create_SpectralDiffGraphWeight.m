function [G, num_edge, info] = Create_SpectralDiffGraphWeight(X, num_segments, k_lap, sigma_l)
% Create_SpectralDiffGraphWeight
%   隣接バンド“差分”ノード (K = N3-1) 上に代表差分の類似度から
%   スペクトルグラフを構築し、重み付きインシデンス行列 G を返す。
%
% 入力:
%   X            : [N1 x N2 x N3] HSI
%   num_segments : セグメント数（imsegkmeans）
%   k_lap        : 相互kNNのk（既定: min(10, K-1)）
%   sigma_l      : "med" → 中央値ベース, "90" → 90%値ベース,
%                  数値 → その値を直接利用
%
% 出力:
%   G        : [K x E] 重み付きインシデンス行列 (gpuArray(single), 疎行列)
%   num_edge : エッジ数
%   info     : 補助情報構造体

[n1, n2, n3] = size(X);
K = n3 - 1;
if nargin < 3 || isempty(k_lap), k_lap = min(10, K-1); end
k_lap = max(0, min(k_lap, K-1));

% --- 1) セグメント化 ---
guide_image = mean(X, 3);
labels = imsegkmeans(guide_image, num_segments);
segID = reshape(labels, n1, n2);

% --- 2) 代表スペクトル ---
S = num_segments;
r = zeros(S, n3);
counts = zeros(S,1);
for s = 1:S
    mask = (segID == s);
    counts(s) = nnz(mask);
    if counts(s) > 0
        tmp = X .* repmat(mask, [1,1,n3]);
        r(s,:) = sum(reshape(tmp, [], n3), 1) / counts(s);
    end
end

% --- 3) 隣接バンド差分 ---
Rdiff = r(:,2:end) - r(:,1:end-1);   % [S x K]

% --- 4) 距離^2行列 ---
Wseg = diag(counts + eps);
Gmat = sqrt(Wseg) * Rdiff; 
GTG  = Gmat.' * Gmat;
nrm2 = diag(GTG);
D2   = (nrm2 + nrm2.') - 2*GTG;
D2(1:K+1:end) = 0;

% --- 5) sigma_l の決定 ---
d2_vec = D2(triu(true(K),1));
d2_vec = d2_vec(d2_vec > 0);

if ischar(sigma_l) || isstring(sigma_l)
    switch lower(string(sigma_l))
        case "med"
            val_sigma_l = sqrt(median(d2_vec));
            mode_sigma  = "median";
        case "90"
            val_sigma_l = sqrt(prctile(d2_vec,90));
            mode_sigma  = "p90";
        otherwise
            error('sigma_l must be "med", "90", or a numeric value.');
    end
elseif isnumeric(sigma_l)
    val_sigma_l = sigma_l;
    mode_sigma  = "manual";
else
    error('sigma_l must be "med", "90", or a numeric value.');
end

% --- 6) 類似度 W_band ---
W_band = exp(-D2./(2*val_sigma_l^2));
W_band(1:K+1:end) = 0;

% --- 7) 相互kNNで疎化 ---
if k_lap > 0 && K > 1
    Wk = zeros(K);
    for i = 1:K
        [~, idx] = maxk(W_band(i,:), min(k_lap, K-1));
        Wk(i, idx) = W_band(i, idx);
    end
    W_band = max(Wk, Wk.');
end

% --- 8) インシデンス行列 G ---
% （数値安定のため W_band を double にしてから find ）
[i_idx, j_idx, wij] = find(triu(double(W_band), 1));   % E edges
num_edge = numel(wij);

if num_edge == 0
    Gcpu = sparse(K, 0);
else
    % 列インデックスを2回並べて 2E にする（各列=1本のエッジに +w と -w を立てる）
    ii = [i_idx; j_idx];                 % 2E x 1
    jj = [(1:num_edge)'; (1:num_edge)']; % 2E x 1
    ss = [wij; -wij];                    % 2E x 1
    Gcpu = sparse(ii, jj, ss, K, num_edge);
end
G = gpuArray(Gcpu);

% --- 9) 情報まとめ ---
info = struct();
info.K              = K;
info.labels         = segID;
info.representative = r;
info.Rdiff          = Rdiff;
info.segment_sizes  = counts;
info.W_band         = W_band;
info.L_delta        = diag(sum(W_band,2)) - W_band;
info.sigma_l_mode   = mode_sigma;
info.val_sigma_l    = val_sigma_l;
info.k_lap          = k_lap;
info.edge_ijw       = [i_idx, j_idx, wij];

end
