%% visualize_param_surface_GASSTV_spagraph.m
% This script loads val_mpsnr and val_mssim saved by
% exp_denoising_gstd_GASSTV_param_spagraph.m (your format),
% and visualizes the relation between (sigma_sp, lambda_rho_sp) and metrics.

clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

%% Selecting conditions (same as your experiment)
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

images = { ...
    "JasperRidge", ...
    "PaviaU", ...
};
idc_images = 1:numel(images);

load("dir_save_comp_folder.mat", "dir_save_comp_folder")


is_output = 1;


%% Setting common parameters
rhos = {0.95};
rho_radius = 0.95;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;

%% Method info (same as your experiment)
ASSTV_CondatVu.sigma_sp = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
GASSTV_CondatVu.lambda_rho_sp = [0.7, 0.8, 0.9, 1, 1.1, 1.2];

GASSTV_CondatVu.sigma_l = {"med"};
GASSTV_CondatVu.lambda2 = [1];
GASSTV_CondatVu.k_lap = [10];
GASSTV_CondatVu.num_segments = [4];
GASSTV_CondatVu.order_filt = [5];

name_method = "GASSTV_CondatVu";

get_params_savetext = @(params) ...
    sprintf("l%.2g_%.2g_sig%g_%s_ns%d_fil%d_k%g_r%.2f_stop1e-%d", ...
    params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
    params.num_segments, params.order_filt, params.k_lap, ...
    rho_radius, stopcri_idx);

%% ===== MAIN LOOP: load metrics into grids and visualize =====
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images

    %% Generating observation (to reproduce folder name + epsilon/alpha/beta if you need later)
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
    image = images{idx_image};

    [HSI_clean, hsi] = Load_HSI(image);
    noise_seed = "default";
    [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

    [~, ~, val_sigma_sp_med] = Create_SpatialGraphWeight(HSI_noisy, "med");
    [~, ~, val_sigma_sp_90] = Create_SpatialGraphWeight(HSI_noisy, "90");


    % Parameter axes you want to analyze
    sigma_list  = GASSTV_CondatVu.sigma_sp(:);
    lambda_list = GASSTV_CondatVu.lambda_rho_sp(:);

    % Result grids (rows=sigma_sp, cols=lambda_rho_sp)
    vals_MPSNR = nan(numel(sigma_list), numel(lambda_list));
    vals_MSSIM = nan(numel(sigma_list), numel(lambda_list));

    %% Loop over only (sigma_sp, lambda_rho_sp)
    % Keep others fixed as in your experiment
    for i = 1:numel(sigma_list)
        for j = 1:numel(lambda_list)

            params = struct();
            params.lambda_rho_sp = lambda_list(j);
            params.lambda2 = GASSTV_CondatVu.lambda2(1);
            params.sigma_sp = sigma_list(i);
            params.sigma_l = GASSTV_CondatVu.sigma_l{1};
            params.num_segments = GASSTV_CondatVu.num_segments(1);
            params.order_filt = GASSTV_CondatVu.order_filt(1);
            params.k_lap = GASSTV_CondatVu.k_lap(1);
            params.maxiter = maxiter;
            params.stopcri = stopcri;
            params.rho_radius = rhos{1};

            % If rho_radius exists, calculate epsilon/alpha/beta (same as your experiment)
            rate_except_sd = 1 - deg.sparse_rate - 2*deg.deadline_rate + 2*deg.sparse_rate*deg.deadline_rate;
            params.epsilon = params.rho_radius * deg.gaussian_sigma * sqrt(hsi.N * rate_except_sd);
            params.alpha   = params.rho_radius * (0.5 * hsi.N * (deg.sparse_rate - 2*deg.sparse_rate*deg.deadline_rate));
            params.beta    = params.rho_radius * hsi.N * rate_except_sd ...
                * deg.stripe_rate * deg.stripe_intensity / 2;

            params_savetext = get_params_savetext(params);

            %% Loading metrics
            dir_result_file = fullfile( ...
                dir_save_comp_folder, ...
                append("denoising_", image), ...
                append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                       "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
                name_method, ...
                append(params_savetext, ".mat") ...
            );

            if ~isfile(dir_result_file)
                % missing result -> leave NaN
                fprintf("[MISS] %s\n", dir_result_file);
                continue;
            end

            S = load(dir_result_file, "val_mpsnr", "val_mssim");

            if isfield(S, "val_mpsnr"); vals_MPSNR(i,j) = S.val_mpsnr; end
            if isfield(S, "val_mssim"); vals_MSSIM(i,j) = S.val_mssim; end
        end
    end

    %% ===== Visualization (heatmaps) =====
    figure('Units','pixels','Position',[100 100 1400 900]);  % 大きめ
    t = tiledlayout(1,2, "TileSpacing","compact", "Padding","compact");
    
    % Axes values
    x = lambda_list(:).';          % omega_sp
    y = log10(sigma_list(:));      % log10(sigma_sp)
    
    % Reference sigma values
    y_med = log10(val_sigma_sp_med);
    y_90  = log10(val_sigma_sp_90);
    
    % -------- Font sizes (≈3–4×) --------
    fs_axis  = 32;   % 目盛り
    fs_label = 37;   % 軸ラベル
    fs_text  = 30;   % 注釈（med / 90）
    
    % ================== MPSNR ==================
    nexttile;
    imagesc(x, y, vals_MPSNR);
    set(gca, ...
        "YDir","normal", ...
        "FontSize", fs_axis, ...
        "XColor","k", ...
        "YColor","k", ...
        "LineWidth", 1.2);
    
    xlabel("\omega_{sp}", "FontSize", fs_label, "Color","k");
    ylabel("log_{10}(\sigma_{sp})", "FontSize", fs_label, "Color","k");
    
    cb1 = colorbar;
    set(cb1, ...
        "Color","k", ...
        "FontSize", fs_axis, ...
        "LineWidth", 1.2);
    
    hold on;
    yline(y_med, "--k", ...
        sprintf(" med : %.3g", val_sigma_sp_med), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    
    yline(y_90, "-.k", ...
        sprintf(" 90%% : %.3g", val_sigma_sp_90), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    hold off;
    
    % ================== MSSIM ==================
    nexttile;
    imagesc(x, y, vals_MSSIM);
    set(gca, ...
        "YDir","normal", ...
        "FontSize", fs_axis, ...
        "XColor","k", ...
        "YColor","k", ...
        "LineWidth", 1.2);
    
    xlabel("\omega_{sp}", "FontSize", fs_label, "Color","k");
    ylabel("log_{10}(\sigma_{sp})", "FontSize", fs_label, "Color","k");
    
    cb2 = colorbar;
    set(cb2, ...
        "Color","k", ...
        "FontSize", fs_axis, ...
        "LineWidth", 1.2);
    
    hold on;
    yline(y_med, "--k", ...
        sprintf(" med : %.3g", val_sigma_sp_med), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    
    yline(y_90, "-.k", ...
        sprintf(" 90%% : %.3g", val_sigma_sp_90), ...
        "LineWidth", 3.5, ...
        "FontSize", fs_text, ...
        "FontWeight","bold");
    hold off;
    
    %% ===== PNG output only =====
    if exist("is_output","var")

        save_figure_folder = fullfile(...
            dir_save_comp_folder, ...
            "GASSTV_Defence", ...
            "ParamAnal");
    
        if ~exist(save_figure_folder,"dir")
            mkdir(save_figure_folder);
        end
    
        fig_name = sprintf( ...
            "%s_GASSTV_g%g_ps%g_pt%g_pd%g_heatmap", ...
            image, ...
            deg.gaussian_sigma, deg.sparse_rate, ...
            deg.stripe_rate, deg.deadline_rate );
    
        fig_path = fullfile(save_figure_folder, fig_name + ".png");
    
        set(gcf,"Color","w");
        print(gcf, fig_path, "-dpng", "-r500");   % 500 dpi（論文向け）
    
        fprintf("[Saved PNG] %s\n", fig_path);
    end


    %% ===== Optional: 3D surfaces =====
    % Uncomment if you want surfaces
    %{
    [X, Y] = meshgrid(x, y);

    figure("Name", sprintf("%s MPSNR surf", image));
    surf(X, Y, MPSNR, "EdgeColor","none");
    xlabel("\lambda_{\rho,sp}"); ylabel("log_{10}(\sigma_{sp})"); zlabel("MPSNR");
    title(sprintf("%s: MPSNR surface", image), "Interpreter","none");
    colorbar; view(45,30);

    figure("Name", sprintf("%s MSSIM surf", image));
    surf(X, Y, MSSIM, "EdgeColor","none");
    xlabel("\lambda_{\rho,sp}"); ylabel("log_{10}(\sigma_{sp})"); zlabel("MSSIM");
    title(sprintf("%s: MSSIM surface", image), "Interpreter","none");
    colorbar; view(45,30);
    %}
end
end

fprintf("******* finis *******\n");
