%% visualize_param_surface_GASSTV_sigmaL_vs_segments.m
% Heatmaps for (num_segments, sigma_l) vs (MPSNR, MSSIM)
% Styling aligned with your spagraph-param heatmaps.

clear
close all; clc;

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
    {0.1,   0.05,  0.05,  0.5,  0},     ... % <- used
    {0.1,   0,     0,     0,    0.01},  ...
    {0.1,   0.05,  0.05,  0.5,  0.01},  ...
};

idc_noise_conditions = [7];

images = { "JasperRidge", "PaviaU" };
idc_images = 1:numel(images);

load("dir_save_comp_folder.mat", "dir_save_comp_folder")

is_output = 1;

%% Common parameters
rhos = {0.95};
rho_radius = 0.95;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;

%% Method info
GASSTV_CondatVu.sigma_l       = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
GASSTV_CondatVu.num_segments  = [2, 3, 4, 5, 6, 7, 8, 9, 10];
% GASSTV_CondatVu.num_segments  = [3,5,7,9];

GASSTV_CondatVu.sigma_sp      = {"med", "90"};   % image-dependent (keep as-is)
GASSTV_CondatVu.lambda_rho_sp = [1];
GASSTV_CondatVu.lambda2       = [1];
GASSTV_CondatVu.k_lap         = [1000];
GASSTV_CondatVu.order_filt    = [5];

name_method = "GASSTV_CondatVu";

% NOTE: sigma_sp is %s, sigma_l is %g (your requirement)
get_params_savetext = @(params) ...
    sprintf("l%.2g_%.2g_sig%s_%g_ns%d_fil%d_r%.2f_stop1e-%d", ...
    params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
    params.num_segments, params.order_filt, ...
    rho_radius, stopcri_idx);

%% ===== MAIN LOOP =====
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images

    %% Observation (for folder naming consistency)
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

    %% Axes lists (FIXED)
    seg_list     = GASSTV_CondatVu.num_segments(:); % x (integer)
    sigma_l_list = GASSTV_CondatVu.sigma_l(:);      % y

    % Z must be [numel(y) x numel(x)] for imagesc(x,y,Z)
    vals_MPSNR = nan(numel(sigma_l_list), numel(seg_list));
    vals_MSSIM = nan(numel(sigma_l_list), numel(seg_list));

    %% Reference sigma_l values (med/90) computed from current noisy HSI
    % (The function returns val_sigma_l; depends on your implementation)
    [~, ~, ~, info_l] = Create_SpectralGraphLaplacian(HSI_noisy, 4, "med", 10, 5);
    val_sigma_l_med = info_l.val_sigma_l;

    [~, ~, ~, info_l] = Create_SpectralGraphLaplacian(HSI_noisy, 4, "90", 10, 5);
    val_sigma_l_90 = info_l.val_sigma_l;

    %% Load metrics
    for iy = 1:numel(sigma_l_list)       % y index
        for ix = 1:numel(seg_list)       % x index

            params = struct();
            params.lambda_rho_sp = GASSTV_CondatVu.lambda_rho_sp(1);
            params.lambda2       = GASSTV_CondatVu.lambda2(1);
            params.sigma_sp      = GASSTV_CondatVu.sigma_sp{idx_image};
            params.sigma_l       = sigma_l_list(iy);
            params.num_segments  = seg_list(ix);
            params.order_filt    = GASSTV_CondatVu.order_filt(1);
            params.k_lap         = GASSTV_CondatVu.k_lap(1);
            params.maxiter       = maxiter;
            params.stopcri       = stopcri;
            params.rho_radius    = rhos{1};

            params_savetext = get_params_savetext(params);

            dir_result_file = fullfile( ...
                dir_save_comp_folder, ...
                append("denoising_", image), ...
                append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                       "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
                name_method, ...
                append(params_savetext, ".mat") ...
            );

            if ~isfile(dir_result_file)
                fprintf("[MISS] %s\n", dir_result_file);
                continue;
            end

            S = load(dir_result_file, "val_mpsnr", "val_mssim");
            if isfield(S,"val_mpsnr"); vals_MPSNR(iy,ix) = S.val_mpsnr; end
            if isfield(S,"val_mssim"); vals_MSSIM(iy,ix) = S.val_mssim; end
        end
    end

    %% ===== Visualization (match spagraph style) =====
    fig = figure("Color","w", "Position",[120 120 1700 650]);
    tiledlayout(1,2, "TileSpacing","compact", "Padding","compact");

    % x is integer segments, y is log10(sigma_l)
    x = seg_list(:).';
    y = log10(sigma_l_list(:));

    y_med = log10(val_sigma_l_med);
    y_90  = log10(val_sigma_l_90);

    % Font sizes (adjust to your caption-sized style)
    fs_axis  = 32;
    fs_label = 30;
    fs_text  = 37;

    % ---------- MPSNR ----------
    nexttile;
    imagesc(x, y, vals_MPSNR);
    set(gca, "YDir","normal", "FontSize",fs_axis, "XColor","k", "YColor","k", ...
        "LineWidth", 1.2);
    xlabel("K", "FontSize",fs_label, "Color","k");             % S = num_segments
    ylabel("log_{10}(\sigma_{l})", "FontSize",fs_label, "Color","k");
    xticks(x); % show all segment counts

    cb1 = colorbar; set(cb1, "Color","k", "FontSize",fs_axis, "Box","on");

    hold on;
    yline(y_med, "--k", ...
        sprintf(" med : %.3g", val_sigma_l_med), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    
    yline(y_90, "-.k", ...
        sprintf(" 90%% : %.3g", val_sigma_l_90), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    hold off;

    % ---------- MSSIM ----------
    nexttile;
    imagesc(x, y, vals_MSSIM);
    set(gca, "YDir","normal", "FontSize",fs_axis, "XColor","k", "YColor","k", ...
        "LineWidth", 1.2);
    xlabel("K", "FontSize",fs_label, "Color","k");
    ylabel("log_{10}(\sigma_{l})", "FontSize",fs_label, "Color","k");
    xticks(x);

    cb2 = colorbar; set(cb2, "Color","k", "FontSize",fs_axis, "Box","on");

    hold on;
    yline(y_med, "--k", ...
        sprintf(" med : %.3g", val_sigma_l_med), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    
    yline(y_90, "-.k", ...
        sprintf(" 90%% : %.3g", val_sigma_l_90), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    hold off;

    %% ===== Output (PNG) =====
    if exist("is_output", "var")
        save_figure_folder = fullfile( ...
            dir_save_comp_folder, ...
            "GASSTV_Defence", ...
            "ParamAnal");

        if ~exist(save_figure_folder, "dir")
            mkdir(save_figure_folder);
        end

        fig_base = sprintf("%s_GASSTV_g%g_ps%g_pt%g_pd%g_sigmaL_vs_segments", ...
            image, deg.gaussian_sigma, deg.sparse_rate, deg.stripe_rate, deg.deadline_rate);

        fig_path_png = fullfile(save_figure_folder, fig_base + ".png");
        print(fig, fig_path_png, "-dpng", "-r200");
        fprintf("[Saved PNG] %s\n", fig_path_png);
    end

end
end

fprintf("******* finis *******\n");
