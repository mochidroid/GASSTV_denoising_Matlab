%% visualize_param_curve_GASSTV_klap_with_edges.m
% - Loads val_mpsnr / val_mssim from saved results
% - Computes edge counts using Create_SpectralGraphLaplacian for each k_lap (=M)
% - Plots: [MPSNR vs M] [MSSIM vs M] [#edges vs M] in one figure
%
% Policy:
% - sigma_sp is selected depending on image index (as you requested)
% - sigma_l is fixed as "med" (string) in this script
%
% Output:
% - PNG only, triggered by existence of variable 'is_output' (comment out to disable)

clear
close all
clc

addpath(genpath("sub_functions"))
addpath("func_metrics")

%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ...
    {0.05,  0.05,  0,     0,    0},     ...
    {0.1,   0.05,  0,     0,    0},     ...
    {0.05,  0,     0.05,  0.5,  0},     ...
    {0.1,   0,     0.05,  0.5,  0},     ...
    {0.05,  0.05,  0.05,  0.5,  0},     ...
    {0.1,   0.05,  0.05,  0.5,  0},     ... % <- you used this
    {0.1,   0,     0,     0,    0.01},  ...
    {0.1,   0.05,  0.05,  0.5,  0.01},  ...
};
idc_noise_conditions = [7];

images = {"JasperRidge", "PaviaU"};
idc_images = 1:numel(images);

load("dir_save_comp_folder.mat", "dir_save_comp_folder")

% is_output = 1;   % If this variable exists, figures are saved (comment out to disable)

%% Common parameters
rhos = {0.95};
rho_radius = 0.95;

stopcri_idx = 5;
stopcri = 10^(-stopcri_idx);

maxiter = 20000;

%% Method info (sweep only k_lap)
GASSTV_CondatVu.sigma_l = {0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1};                     % fixed (string)
GASSTV_CondatVu.k_lap   = [2, 5, 8, 10, 15, 20, 30, 40];  % sweep

GASSTV_CondatVu.sigma_sp = {"med", "90"};              % image-dependent selection
GASSTV_CondatVu.lambda_rho_sp = [0.9];                 % fixed
GASSTV_CondatVu.lambda2 = [1];                         % fixed
GASSTV_CondatVu.num_segments = [4];                    % fixed
GASSTV_CondatVu.order_filt = [5];                      % fixed

name_method = "GASSTV_CondatVu";

% IMPORTANT: match the experiment script's naming exactly.
% sigma_sp and sigma_l are strings -> use %s for both.
get_params_savetext = @(params) ...
    sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_k%g_r%.2f_stop1e-%d", ...
    params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
    params.num_segments, params.order_filt, params.k_lap, ...
    rho_radius, stopcri_idx);

%% ===== MAIN LOOP =====
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images

    %% Observation info (folder naming + edge computation needs HSI_noisy)
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
    image = images{idx_image};

    [HSI_clean, hsi] = Load_HSI(image);
    noise_seed = "default";
    [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

    HSI_noisy = single(HSI_noisy);
    clear HSI_clean hsi

    %% Sweep axis
    M_list = GASSTV_CondatVu.k_lap(:).';   % k_lap sweep

    % Metrics
    vals_MPSNR = nan(size(M_list));
    vals_MSSIM = nan(size(M_list));

    % Edge counts over segments
    E_mean = nan(size(M_list));
    E_min  = nan(size(M_list));
    E_max  = nan(size(M_list));

    % Tags
    sigsp_tag = string(GASSTV_CondatVu.sigma_sp{idx_image}); % image-dependent
    sigl_tag  = string(GASSTV_CondatVu.sigma_l{1});          % fixed "med"

    num_segments = GASSTV_CondatVu.num_segments(1);
    order_filt   = GASSTV_CondatVu.order_filt(1);

    %% Loop over M: load metrics + compute edges
    for j = 1:numel(M_list)

        M = M_list(j);

        % ---- build params (for file name) ----
        params = struct();
        params.lambda_rho_sp = GASSTV_CondatVu.lambda_rho_sp(1);
        params.lambda2       = GASSTV_CondatVu.lambda2(1);
        params.sigma_sp      = GASSTV_CondatVu.sigma_sp{idx_image};  % image-dependent
        params.sigma_l       = GASSTV_CondatVu.sigma_l{1};           % string
        params.num_segments  = num_segments;
        params.order_filt    = order_filt;
        params.k_lap         = M;
        params.maxiter       = maxiter;
        params.stopcri       = stopcri;
        params.rho_radius    = rhos{1};

        params_savetext = get_params_savetext(params);

        dir_result_file = fullfile( ...
            dir_save_comp_folder, ...
            "denoising_" + image, ...
            "g" + num2str(deg.gaussian_sigma) + "_ps" + num2str(deg.sparse_rate) + ...
            "_pt" + num2str(deg.stripe_rate) + "_pd" + num2str(deg.deadline_rate), ...
            name_method, ...
            params_savetext + ".mat" ...
        );

        % ---- load metrics ----
        if isfile(dir_result_file)
            S = load(dir_result_file, "val_mpsnr", "val_mssim");
            if isfield(S,"val_mpsnr"); vals_MPSNR(j) = S.val_mpsnr; end
            if isfield(S,"val_mssim"); vals_MSSIM(j) = S.val_mssim; end
        else
            fprintf("[MISS] %s\n", dir_result_file);
        end

        % ---- compute edge counts (undirected) ----
        % Uses W_band_all(:,:,s) returned in info_l
        try
            [~, ~, ~, info_l] = Create_SpectralGraphLaplacian( ...
                HSI_noisy, num_segments, sigl_tag, M, order_filt);

            W_all = info_l.W_band_all;   % [K x K x S]
            Sseg = size(W_all, 3);

            E_seg = zeros(Sseg,1);
            for s = 1:Sseg
                Ws = W_all(:,:,s);
                E_seg(s) = nnz(triu(Ws > 0, 1));  % count i<j
            end

            E_mean(j) = mean(E_seg);
            E_min(j)  = min(E_seg);
            E_max(j)  = max(E_seg);
        catch ME
            fprintf("[EDGE-FAIL] M=%g : %s\n", M, ME.message);
        end
    end

    %% ===== Integrated visualization (3 panels) =====
    figure('Units','pixels','Position',[100 100 2100 650]);
    t = tiledlayout(1,3, "TileSpacing","compact", "Padding","compact");

    % Large fonts (caption-level)
    fs_axis  = 30;
    fs_label = 34;

    % --- Panel 1: MPSNR vs M ---
    nexttile;
    plot(M_list, vals_MPSNR, "-o", "LineWidth", 2.5, "MarkerSize", 8);
    grid on;
    set(gca, "FontSize", fs_axis, "XColor","k","YColor","k","LineWidth",1.2);
    xlabel("M (k_{lap})", "FontSize", fs_label, "Color","k");
    ylabel("MPSNR", "FontSize", fs_label, "Color","k");

    % --- Panel 2: MSSIM vs M ---
    nexttile;
    plot(M_list, vals_MSSIM, "-o", "LineWidth", 2.5, "MarkerSize", 8);
    grid on;
    set(gca, "FontSize", fs_axis, "XColor","k","YColor","k","LineWidth",1.2);
    xlabel("M (k_{lap})", "FontSize", fs_label, "Color","k");
    ylabel("MSSIM", "FontSize", fs_label, "Color","k");

    % --- Panel 3: #edges vs M (mean + min-max over segments) ---
    nexttile;
    % mean
    plot(M_list, E_mean, "-o", "LineWidth", 2.5, "MarkerSize", 8); hold on;

    % min-max bars
    for j = 1:numel(M_list)
        if ~isnan(E_min(j)) && ~isnan(E_max(j))
            plot([M_list(j) M_list(j)], [E_min(j) E_max(j)], "-", "LineWidth", 2.0);
        end
    end
    hold off;

    grid on;
    set(gca, "FontSize", fs_axis, "XColor","k","YColor","k","LineWidth",1.2);
    xlabel("M (k_{lap})", "FontSize", fs_label, "Color","k");
    ylabel("#edges (undirected)", "FontSize", fs_label, "Color","k");

    %% ===== Output PNG only =====
    if exist("is_output","var")
        save_figure_folder = fullfile(dir_save_comp_folder, "GASSTV_Defence", "ParamAnal");
        if ~exist(save_figure_folder, "dir")
            mkdir(save_figure_folder);
        end

        fig_base = sprintf("%s_%s_sigsp%s_sigl%s_g%g_ps%g_pt%g_pd%g_M_metrics_edges", ...
            image, name_method, sigsp_tag, sigl_tag, ...
            deg.gaussian_sigma, deg.sparse_rate, deg.stripe_rate, deg.deadline_rate);

        fig_path_png = fullfile(save_figure_folder, fig_base + ".png");

        set(gcf, "Color","w");
        print(gcf, fig_path_png, "-dpng", "-r500");
        fprintf("[Saved PNG] %s\n", fig_path_png);
    end

end
end

fprintf("******* finis *******\n");
