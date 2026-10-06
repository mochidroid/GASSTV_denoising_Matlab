function [L_delta, B, K, info_l] = Create_GraphWeight_GLR(X, opts)
% Create_SpectralDiffGraphLaplacian
% -------------------------------------------------------------------------
% Build a *spectral-difference* graph over K = N3-1 nodes (each node = an
% adjacent-band difference) using segment-wise representative spectra, then
% return:
%   - L_delta           : [K x K] graph Laplacian (gpuArray-single, symmetric)
%   - B                 : [K x K] matrix s.t. B = L_delta^{1/2} (gpuArray-single)
%   - lam_max_Ldelta    : largest eigenvalue of L_delta (single)
%   - info_l            : struct with construction details (including CPU copies)
%
% INPUTS:
%   X    : [N1 x N2 x N3] HSI (guide or observed)
%   opts : struct with fields
%       .num_segments  (default = 50)   % #segments for imsegkmeans
%       .k_lap         (default = min(10, K-1)) % mutual kNN sparsification
%       .sigma_l       (default = "med") % "med" | "90" | numeric (Gaussian scale)
%
% NOTES:
%   - Distances between nodes are computed from *segment-weighted* distances
%     of representative *adjacent-band differences*.
%   - We build an affinity W_band via Gaussian kernel and mutual kNN,
%     then L = D - W.
%   - B and lam_max are computed on CPU (double) for robustness, then moved to GPU.
% -------------------------------------------------------------------------

[n1, n2, n3] = size(X);
K = n3 - 1;

% ---- defaults ----
if ~isfield(opts,'num_segments'), opts.num_segments = 50; end
if ~isfield(opts,'k_lap'),        opts.k_lap        = min(10, max(0,K-1)); end
if ~isfield(opts,'sigma_l'),      opts.sigma_l      = "med"; end

S      = opts.num_segments;
k_lap  = max(0, min(opts.k_lap, max(0,K-1)));
sigma_l_opt = opts.sigma_l;

% ---- 1) Segmentation on a guide image (mean over bands) ----
guide = mean(X,3);
labels = imsegkmeans(guide, S);
segID  = reshape(labels, n1, n2);

% ---- 2) Segment representative spectra r(s,:) ----
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

% ---- 3) Adjacent-band differences Rdiff(s, :) ----
Rdiff = r(:,2:end) - r(:,1:end-1);  % [S x K]

% ---- 4) Weighted squared-distance matrix D2 (K x K) ----
Wseg = diag(counts + eps);           % segment weights = #pixels
Gmat = sqrt(Wseg) * Rdiff;           % [S x K]
GTG  = Gmat.' * Gmat;                % [K x K]
nrm2 = diag(GTG);
D2   = (nrm2 + nrm2.') - 2*GTG;      % pairwise squared distances
D2(1:K+1:end) = 0;

% ---- 5) Choose sigma_l: "med" / "90" / numeric ----
d2_vec = D2(triu(true(K),1));
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

% ---- 6) Affinity W_band and mutual kNN sparsification ----
W_band = exp( -double(D2) / (2*val_sigma_l^2) );
W_band(1:K+1:end) = 0;

if k_lap > 0 && K > 1
    Wk = zeros(K);
    for i = 1:K
        [~, idx] = maxk(W_band(i,:), min(k_lap, K-1));
        Wk(i, idx) = W_band(i, idx);
    end
    W_band = max(Wk, W_band.');  % symmetric mutual-kNN
end

% ---- 7) Laplacian L = D - W ----
d = sum(W_band, 2);
L_cpu = diag(d) - W_band;      % double on CPU for eig robustness

% ---- 8) B = L^{1/2} (CPU eig) and lam_max(L) ----
% Symmetrize numerically
L_sym = (L_cpu + L_cpu.')/2;
% Eigen-decomposition
[V, D2_eig] = eig(L_sym);
evals = real(diag(D2_eig));
evals(evals < 0) = 0;                 % clip tiny negatives
B_cpu  = V * diag(sqrt(evals)) * V.'; % L^{1/2}
% Largest eigenvalue
try
    lam_max_Ldelta = eigs(L_sym, 1, 'lm');   % double
    lam_max_Ldelta = real(lam_max_Ldelta);
catch
    lam_max_Ldelta = max(real(eig(L_sym)));
end
lam_max_Ldelta = single(max(lam_max_Ldelta, 0));

% ---- 9) Move to GPU (single) ----
L_delta = gpuArray(single(L_cpu));
B       = gpuArray(single(B_cpu));

% ---- 10) info_l ----
% also provide an edge list (upper-tri) for debug/plot if needed
[i_idx, j_idx, wij] = find(triu(W_band, 1));
info_l = struct();
info_l.K                 = K;
info_l.num_segments      = S;
info_l.segment_sizes     = counts;
info_l.labels            = segID;
info_l.representative    = r;         % [S x N3]
info_l.Rdiff             = Rdiff;     % [S x K]
info_l.W_band            = W_band;    % [K x K]
info_l.sigma_l_mode      = sigma_l_mode;
info_l.val_sigma_l       = val_sigma_l;
info_l.k_lap             = k_lap;
info_l.edge_ijw          = [i_idx, j_idx, wij];   % undirected edges
% CPU copies (for diagnostics/reuse)
info_l.L_delta_cpu       = L_cpu;     % double
info_l.B_cpu             = B_cpu;     % double
info_l.lam_max_Ldelta    = lam_max_Ldelta;  % single (also returned)
end
