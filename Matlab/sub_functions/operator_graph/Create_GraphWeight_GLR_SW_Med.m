function [L_delta, B, lam_max_Ldelta, info_l] = Create_GraphWeight_GLR_SW_Med(X, opts)
% Create_GraphWeight_GLR_SW_Med
% segment-wise (representative-spectrum-wise) spectral GLR graph
%
% INPUT:
%   X    : [N1 x N2 x N3] HSI (guide or clean)
%   opts :
%       .num_segments  (default = 50)
%       .k_lap         (default = min(10, K-1))
%       .sigma_l       (default = "med")  % "med" | "90" | numeric
%
% OUTPUT:
%   L_delta       : [K x K x S] Laplacians (double, CPU)
%   B             : [K x K x S] s.t. B(:,:,s) = L_delta(:,:,s)^{1/2}
%   lam_max_Ldelta: [S x 1] 最大固有値 (single)
%   info_l        : 構築情報（segID, representative spectra など）

[n1, n2, n3] = size(X);
K = n3 - 1;

% ---- defaults ----
if ~isfield(opts,'num_segments'), opts.num_segments = 5; end
if ~isfield(opts,'k_lap'),        opts.k_lap        = min(10, max(0,K-1)); end
if ~isfield(opts,'sigma_l'),      opts.sigma_l      = "med"; end
if ~isfield(opts, 'order_filt'), opts.order_filt = 3; end

S      = opts.num_segments;
k_lap  = max(0, min(opts.k_lap, max(0,K-1)));
sigma_l_opt = opts.sigma_l;
order_filt = opts.order_filt;

% ---- 1) Segmentation on a guide image (mean over bands) ----
guide  = mean(X,3);
labels = imsegkmeans(guide, S);
segID  = reshape(labels, n1, n2);

% ---- 2) Segment representative spectra r(s,:) ----
r      = zeros(S, n3);
counts = zeros(S,1);
for s = 1:S
    mask = (segID == s);
    counts(s) = nnz(mask);
    if counts(s) > 0
        tmp = X .* repmat(mask, [1,1,n3]);
        r(s,:) = sum(reshape(tmp, [], n3), 1) / counts(s);
        r(s,:) = medfilt1(r(s,:), order_filt);
    end
end



% ---- 3) Adjacent-band differences Rdiff(s, :) ----
Rdiff = r(:,2:end) - r(:,1:end-1);  % [S x K]

% ---- 4) Global D2 from all segments (for common sigma_l) ----
Wseg = diag(counts + eps);           % segment weights
Gmat = sqrt(Wseg) * Rdiff;           % [S x K]
GTG  = Gmat.' * Gmat;                % [K x K]
nrm2 = diag(GTG);
D2_global   = (nrm2 + nrm2.') - 2*GTG;
D2_global(1:K+1:end) = 0;

% ---- 5) Choose sigma_l from GLOBAL distances ----
d2_vec = D2_global(triu(true(K),1));
d2_vec = d2_vec(d2_vec > 0);
if isempty(d2_vec), d2_vec = 1; end

if ischar(sigma_l_opt) || isstring(sigma_l_opt)
    switch lower(string(sigma_l_opt))
        case "med"
            val_sigma_l = sqrt(median(double(d2_vec)));
            sigma_l_mode = "median";
        case "90"
            val_sigma_l = sqrt(prctile(double(d2_vec), 90));
            sigma_l_mode = "p90";
        otherwise
            error('opts.sigma_l must be "med", "90", or a numeric value.');
    end
else
    val_sigma_l  = double(sigma_l_opt);
    sigma_l_mode = "manual";
end
val_sigma_l = max(val_sigma_l, eps);

% ---- 6) per-segment GLR graph ----
L_cpu_all    = zeros(K, K, S);   % double
B_cpu_all    = zeros(K, K, S);   % double
lam_max_all  = zeros(S, 1, 'single');
W_band_all   = zeros(K, K, S);   % for debug/info

for s = 1:S
    diff_s = Rdiff(s,:);          % [1 x K]
    
    % (a) segment-wise squared distances (K x K)
    v = diff_s(:);
    nrm2_s = v.^2;
    D2_s   = (nrm2_s + nrm2_s.') - 2*(v * v.');
    D2_s(1:K+1:end) = 0;
    
    % (b) Gaussian affinities with GLOBAL sigma_l
    W_s = exp( -double(D2_s) / (2*val_sigma_l^2) );
    W_s(1:K+1:end) = 0;
    
    % (c) mutual kNN sparsification
    if k_lap > 0 && K > 1
        Wk = zeros(K);
        for i = 1:K
            [~, idx] = maxk(W_s(i,:), min(k_lap, K-1));
            Wk(i, idx) = W_s(i, idx);
        end
        W_s = max(Wk, W_s.');  % symmetric
    end
    
    % (d) Laplacian L_s = D - W_s
    d = sum(W_s, 2);
    L_s = diag(d) - W_s;
    L_s = (L_s + L_s.')/2;      % symmetrize
    
    % (e) B_s = L_s^{1/2}
    [V, D2_eig] = eig(L_s);
    evals = real(diag(D2_eig));
    evals(evals < 0) = 0;
    B_s  = V * diag(sqrt(evals)) * V.';
    
    % % (f) max eigenvalue
    % try
    %     lam_max = eigs(L_s, 1, 'lm');
    %     lam_max = real(lam_max);
    % catch
    %     lam_max = max(real(eig(L_s)));
    % end
    % lam_max = single(max(lam_max, 0));

    % --- (f) max eigenvalue (robust, no eigs) ---
    eigvals = eig(L_s);                   % 全固有値（複素でもOK）
    eigvals = real(eigvals);              % 実部だけを使う
    eigvals(~isfinite(eigvals)) = 0;      % NaN, Inf を 0 に
    lam_max = max(eigvals);               % 最大値
    
    if ~isfinite(lam_max) || lam_max < 0
        lam_max = 0;                      % 念のためクリップ
    end
    lam_max_all(s) = single(lam_max);
    
    % (g) store
    L_cpu_all(:,:,s)   = L_s;
    B_cpu_all(:,:,s)   = B_s;
    lam_max_all(s)     = lam_max;
    W_band_all(:,:,s)  = W_s;
end

L_delta        = L_cpu_all;      % [K x K x S], double
B              = B_cpu_all;      % [K x K x S], double
lam_max_Ldelta = lam_max_all;    % [S x 1], single

% info for debug
info_l = struct();
info_l.K              = K;
info_l.num_segments   = S;
info_l.order_filt     = order_filt;
info_l.segment_sizes  = counts;
info_l.labels         = segID;          % [n1 x n2]
info_l.representative = r;              % [S x N3]
info_l.Rdiff          = Rdiff;          % [S x K]
info_l.W_band_all     = W_band_all;     % [K x K x S]
info_l.sigma_l_mode   = sigma_l_mode;
info_l.val_sigma_l    = val_sigma_l;
info_l.k_lap          = k_lap;
info_l.L_delta_cpu    = L_cpu_all;
info_l.B_cpu          = B_cpu_all;
end
