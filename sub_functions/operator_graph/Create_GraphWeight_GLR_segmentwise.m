function [L_cells, B_cells, lam_max_vec, info_l] = Create_GraphWeight_GLR_segmentwise(X, opts)
% Create_GraphWeight_GLR_segmentwise
% -------------------------------------------------------------------------
% Build a *segment-wise* spectral-difference graph for each segment s:
%   - Nodes: K = N3 - 1 (adjacent-band differences)
%   - For each segment s, construct a K x K Laplacian L_s from the
%     representative adjacent-band differences of that segment only.
%
% OUTPUT:
%   L_cells    : 1 x S cell, L_cells{s} is [K x K] Laplacian (double)
%   B_cells    : 1 x S cell, B_cells{s} is [K x K] s.t. B_s = L_s^{1/2} (double)
%   lam_max_vec: [S x 1] largest eigenvalue of each L_s (single)
%   info_l     : struct with construction details
%
% INPUT:
%   X    : [n1 x n2 x n3] HSI
%   opts : struct
%       .num_segments  (default = 50)          % #segments for imsegkmeans
%       .k_lap         (default = min(10,K-1)) % mutual kNN sparsification
%       .sigma_l       (default = "med")       % "med" | "90" | numeric
%
% NOTES:
%   - We first compute representative spectra r_s per segment.
%   - For each segment s, we form the 1 x K vector of adjacent-band
%     differences Rdiff(s,:), then build a K x K distance matrix:
%         D2_s(i,j) = (Rdiff_s(i) - Rdiff_s(j))^2
%   - sigma_l is chosen *globally* over all segments (using all distances),
%     then each segment's graph is built with that common scale.
% -------------------------------------------------------------------------

[n1, n2, n3] = size(X);
K = n3 - 1;                    % #nodes = #adjacent-band differences

% ---- defaults ----
if ~isfield(opts,'num_segments'), opts.num_segments = 50; end
if ~isfield(opts,'k_lap'),        opts.k_lap        = min(10, max(0,K-1)); end
if ~isfield(opts,'sigma_l'),      opts.sigma_l      = "med"; end

if numel(opts.k_lap) > 1
    opts.k_lap = opts.k_lap(1);   % 念のため先頭だけ使う
end

S      = opts.num_segments;
k_lap  = max(0, min(opts.k_lap, max(0,K-1)));
sigma_l_opt = opts.sigma_l;

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
    end
end

% ---- 3) Adjacent-band differences Rdiff(s,:) ----
Rdiff = r(:,2:end) - r(:,1:end-1);   % [S x K]

% ---- 4) Collect all squared distances over all segments for sigma_l ----
d2_all = [];
for s = 1:S
    if counts(s) == 0
        continue; % no pixels in this segment
    end
    v  = Rdiff(s,:);              % 1 x K
    vv = v(:);                    % K x 1
    % pairwise squared distances between nodes i,j:
    % D2_s(i,j) = (v(i) - v(j))^2
    D2_s = (vv.^2 + (vv.^2).') - 2*(vv*vv.');
    D2_s(1:K+1:end) = 0;
    d2_vec_s = D2_s(triu(true(K),1));
    d2_vec_s = d2_vec_s(d2_vec_s > 0);
    if ~isempty(d2_vec_s)
        d2_all = [d2_all; d2_vec_s]; %#ok<AGROW>
    end
end
if isempty(d2_all)
    d2_all = 1;
end

% ---- 5) Choose sigma_l: "med" / "90" / numeric ----
if ischar(sigma_l_opt) || isstring(sigma_l_opt)
    switch lower(string(sigma_l_opt))
        case "med"
            val_sigma_l = sqrt(median(double(d2_all)));
            sigma_l_mode = "median";
        case "90"
            val_sigma_l = sqrt(prctile(double(d2_all), 90));
            sigma_l_mode = "p90";
        otherwise
            error('opts.sigma_l must be "med", "90", or a numeric value.');
    end
else
    val_sigma_l  = double(sigma_l_opt);
    sigma_l_mode = "manual";
end
val_sigma_l = max(val_sigma_l, eps);

% ---- 6) For each segment, build W_s, L_s, B_s ----
L_cells     = cell(1,S);
B_cells     = cell(1,S);
lam_max_vec = zeros(S,1,'single');

edge_ijw_cell = cell(1,S);  % edge list [i j w] per segment

for s = 1:S
    if counts(s) == 0
        % empty segment → zero Laplacian
        L_s = zeros(K); 
        B_s = zeros(K); 
        lam_max_vec(s) = 0;
        L_cells{s}     = zeros(K);
        B_cells{s}     = zeros(K);
        edge_ijw_cell{s} = zeros(0,3);
        continue;
    end

    v  = Rdiff(s,:);          % 1 x K
    vv = v(:);                % K x 1

    % squared-distance matrix for this segment
    D2_s = (vv.^2 + (vv.^2).') - 2*(vv*vv.');
    D2_s(1:K+1:end) = 0;

    % Gaussian affinity
    W_s = exp( -double(D2_s) / (2*val_sigma_l^2) );
    W_s(1:K+1:end) = 0;

    % mutual kNN sparsification
    if k_lap > 0 && K > 1
        Wk = zeros(K);
        for i = 1:K
            [~, idx] = maxk(W_s(i,:), min(k_lap, K-1));
            Wk(i, idx) = W_s(i, idx);
        end
        W_s = max(Wk, W_s.');   % symmetric
    end

    % Laplacian L_s = D - W_s
    d_s  = sum(W_s, 2);
    L_s  = diag(d_s) - W_s;
    L_sym = (L_s + L_s.')/2;

    % eig for B_s = L_s^{1/2}
    [V_s, D_s] = eig(L_sym);
    evals_s = real(diag(D_s));
    evals_s(evals_s < 0) = 0;
    B_s = V_s * diag(sqrt(evals_s)) * V_s.';

    % largest eigenvalue
    try
        lam_max = eigs(L_sym, 1, 'lm');
        lam_max = real(lam_max);
    catch
        lam_max = max(real(eig(L_sym)));
    end
    lam_max_vec(s) = single(max(lam_max, 0));

    L_cells{s} = L_s;
    B_cells{s} = B_s;

    % store edge list for this segment (upper triangle)
    [i_idx, j_idx, wij] = find(triu(W_s,1));
    edge_ijw_cell{s} = [i_idx, j_idx, wij];
end

% ---- 7) info_l summary ----
info_l = struct();
info_l.K              = K;
info_l.num_segments   = S;
info_l.segment_sizes  = counts;
info_l.labels         = segID;
info_l.representative = r;        % [S x n3]
info_l.Rdiff          = Rdiff;    % [S x K]
info_l.sigma_l_mode   = sigma_l_mode;
info_l.val_sigma_l    = val_sigma_l;
info_l.k_lap          = k_lap;
info_l.edge_ijw_cell  = edge_ijw_cell;  % per-segment edges [i j w]
info_l.L_cells        = L_cells;       % CPU double
info_l.B_cells        = B_cells;       % CPU double
info_l.lam_max_vec    = lam_max_vec;
end
