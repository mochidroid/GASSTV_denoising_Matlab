function [G_sp, num_edge, info] = Create_NLSpatialGraphTV(X, opts)
% Create_SpatialNLGraphTV
% -------------------------------------------------------------------------
% Builds a nonlocal spatial graph (shared across bands) from patch
% similarity and returns the weighted incidence matrix G_sp ∈ R^{N×E}
% for Spatial Graph Total Variation (GTV).
%
% Candidate edges:
%   Pixels within a square search window of radius Rs (L∞-norm) around each
%   center pixel. (Upper distance bound only; no lower bound here.)
%
% Edge selection priority:
%   Top-K by smallest patch L2 distance (d²). Since
%       w = exp(-d² / (2σ²))
%   is strictly decreasing in d², this is equivalent to Top-K by *largest
%   similarity w*, provided σ is global.
%
% INPUTS:
%   X    : [n1 × n2 × n3] hyperspectral image (used as guide)
%   opts : struct with fields
%       .patch_rad     (default = 1)    % patch radius r ⇒ D = (2r+1)^2
%       .search_rad    (default = 7)    % half-size of local search window
%       .k_nn          (default = 10)   % #neighbors per pixel (mutual kNN)
%       .sigma_sp      (default = "90") % "med" | "90" | numeric
%       .max_per_pixel (default = [])   % optional cap per pixel
%
% OUTPUTS:
%   G_sp     : [N × E] weighted incidence matrix on GPU (gpuArray-sparse-double)
%              each column e places +w at node i and -w at node j
%   num_edge : number of edges E
%   info     : struct with construction details:
%              .sigma_mode, .val_sigma, .deg2_per_node (for ||G||² upper bound)
%              .ijw = [i j w] edge list (1-based indices)
%              .N = [n1 n2], .k_nn, .patch_rad, .search_rad
% -------------------------------------------------------------------------

[n1,n2,~] = size(X);
N = n1*n2;

% ---- defaults ----
if ~isfield(opts,'patch_rad'),      opts.patch_rad      = 1;   end
if ~isfield(opts,'search_rad'),     opts.search_rad     = 7;   end
if ~isfield(opts,'k_nn'),           opts.k_nn           = 10;  end
if ~isfield(opts,'sigma_sp'),       opts.sigma_sp       = "90";end
if ~isfield(opts,'max_per_pixel'),  opts.max_per_pixel  = [];  end

r  = opts.patch_rad;
Rs = opts.search_rad;
K  = opts.k_nn;

% ---- patch features on the guidance (mean over bands here) ----
guid = mean(single(X),3);
pad  = padarray(guid,[r r],'replicate','both');

% vectorized patches: N × D, D = (2r+1)^2
D = (2*r+1)^2;
patches = zeros(N, D, 'single');
col = 1;
for dy = -r:r
    for dx = -r:r
        blk = pad((1:n1)+r+dy, (1:n2)+r+dx);
        patches(:,col) = blk(:);
        col = col+1;
    end
end

% ---- search inside the window and collect K-NN (by smallest d²) ----
ii = []; jj = []; ww = [];      % we temporarily store d² in ww
all_d2 = [];                    % to estimate σ later

% helper: 1D index -> (y,x)
toYX = @(idx) deal( mod(idx-1,n1)+1, floor((idx-1)/n1)+1 );

for p = 1:N
    [py,px] = toYX(p);

    y1 = max(1, py-Rs); y2 = min(n1, py+Rs);
    x1 = max(1, px-Rs); x2 = min(n2, px+Rs);

    cand = reshape( sub2ind([n1,n2], ...
            repmat((y1:y2)', 1, x2-x1+1), ...
            repmat((x1:x2),  y2-y1+1, 1) ), [], 1);

    cand(cand==p) = [];
    if isempty(cand), continue; end

    Pp = patches(p,:);          % 1 × D
    Pc = patches(cand,:);       % C × D

    % patch-wise squared L2 distances
    dif = bsxfun(@minus, Pc, Pp);
    d2  = sum(dif.^2, 2);

    % cache distances to determine σ globally
    all_d2 = [all_d2; d2]; %#ok<AGROW>

    % choose up to K nearest by d² (≡ top-K by similarity for fixed σ)
    k_use = min(K, numel(cand));
    if ~isempty(opts.max_per_pixel)
        k_use = min(k_use, opts.max_per_pixel);
    end
    if k_use==0, continue; end

    [~,ord] = mink(d2, k_use);
    nbr  = cand(ord);
    d2_k = d2(ord);

    ii = [ii; repmat(p,k_use,1)]; %#ok<AGROW>
    jj = [jj; nbr];               %#ok<AGROW>
    ww = [ww; d2_k];              %#ok<AGROW> % keep d²; weights computed later
end

% ---- set σ_sp: "med" / "90" / numeric ----
d2_use = double(all_d2);
if isempty(d2_use), d2_use = 1; end

if ischar(opts.sigma_sp) || isstring(opts.sigma_sp)
    switch lower(string(opts.sigma_sp))
        case "med"
            val_sigma  = sqrt(median(d2_use));
            mode_sigma = "median";
        case "90"
            val_sigma  = sqrt(prctile(d2_use,90));
            mode_sigma = "p90";
        otherwise
            error('sigma_sp must be "med", "90", or a numeric value.');
    end
else
    val_sigma  = double(opts.sigma_sp);
    mode_sigma = "manual";
end
if val_sigma <= 0, val_sigma = 1; end

% ---- convert d² to similarity weights: w = exp(-d² / (2σ²)) ----
wij = exp( -double(ww) / (2*val_sigma^2) );

% ---- mutual k-NN symmetrization (take max for duplicated pairs) ----
I = [ii; jj]; J = [jj; ii]; W = [wij; wij];
[IJ,~,ic] = unique([I J], 'rows');   % unique undirected edges
Wmax = accumarray(ic, W, [], @max);
i_idx = IJ(:,1); j_idx = IJ(:,2); wij = Wmax;

% remove self-loops
mask  = (i_idx ~= j_idx);
i_idx = i_idx(mask); j_idx = j_idx(mask); wij = wij(mask);

% ---- number of edges ----
num_edge = numel(wij);

% ---- weighted incidence G_sp: column e has +w at i, -w at j ----
if num_edge==0
    Gcpu = sparse(N,0);
else
    e    = (1:num_edge).';
    Gcpu = sparse([i_idx; j_idx], [e; e], [wij; -wij], N, num_edge);
end
% IMPORTANT: GPU sparse is double-only in MATLAB. Do NOT cast to single.
G_sp = gpuArray(Gcpu);

% ---- safe upper bound for ||G||² (for step-size tuning) ----
% ||G||² ≤ 2 * max_i ∑_{e incident to i} w_e²
deg2 = full(sum(Gcpu.^2,2));
normG_sq_upper = 2 * max([deg2; 0]);

% ---- info ----
info = struct();
info.N               = [n1 n2];
info.k_nn            = K;
info.patch_rad       = r;
info.search_rad      = Rs;        % upper distance bound (window half-size)
info.sigma_mode      = mode_sigma;
info.val_sigma       = val_sigma;
info.deg2_per_node   = deg2;
info.normG_sq_upper  = normG_sq_upper;
info.ijw             = [i_idx, j_idx, wij];
end
