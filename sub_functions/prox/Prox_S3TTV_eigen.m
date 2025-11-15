function[z] = Prox_S3TTV_eigen(z, gamma)
% z = gather(z);
[v, h, c, d, b1, b2] = size(z);
nb1 = v / b1; nb2 = h / b2;

for i = 1:nb1
    idc_i = (1:b1) + b1*(i-1);
    for j = 1:nb2                
        idc_j = (1:b2) + b2*(j-1);
        for k = 1:b1
            for l = 1:b2
                M = reshape(z(idc_i, idc_j, :, :, k, l), [b1*b2, c*d]);
                % % SVD soft-threshold
                % [U_, S_, V_] = svd(M, 'econ'); 
                % Sthre = diag(max(0, diag(S_) - gamma));
                % z(idc_i, idc_j, :, :, k, l) ...
                %     = reshape(U_ * Sthre * V_', [b1, b2, c, d]);
                X = prox_nuc_thresh(M, gamma);
                z(idc_i, idc_j, :, :, k, l) = reshape(X, [b1, b2, c, d]);
            end
        end
    end
end

function X = prox_nuc_thresh(M, gamma)
% 厳密 prox を、σ>γ の成分だけで再構成（m << n 向け）

[m,n] = size(M);
% 1) Gram
G = M * M.';                 % m x m
G = (G+G.')/2;               % 対称化（数値安定）

% 2) 固有分解（昇順→降順に並べ替え）
[U, D] = eig(G);
lam = max(0, diag(D));
[sig, idx] = sort(sqrt(lam), 'descend');
U = U(:, idx);

% 3) 有効ランク選定
r = nnz(sig > gamma);
if r == 0
    X = zeros(m,n, 'like', M);
    return;
end
Ur   = U(:, 1:r);
s    = sig(1:r);                     % 長さ r
Sinv = diag(1 ./ s);

% 4) 右特異ベクトル（必要分だけ）
Vr = (M.' * Ur) * Sinv;             % n x r

% 5) ソフト閾値して再構成
s_shrunk = s - gamma;               % >0 保証済み
X = Ur * (diag(s_shrunk) * Vr.');

% （数値安定の微調整：σ≈γ の境界は tol を設けると良い）