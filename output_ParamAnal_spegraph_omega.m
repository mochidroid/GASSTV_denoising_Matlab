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


% is_output = 1;


%% Setting common parameters
rhos = {0.95};
rho_radius = 0.95;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;

%% Method info (same as your experiment)
GASSTV_CondatVu.sigma_l = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
GASSTV_CondatVu.lambda2 = [0.01, 0.05, 0.1, 0.5, 1, 5, 10];


GASSTV_CondatVu.sigma_sp = {"med", "90"};
GASSTV_CondatVu.lambda_rho_sp = [0.9];
% GASSTV_CondatVu.lambda2 = [1];
GASSTV_CondatVu.k_lap = [10];
GASSTV_CondatVu.num_segments = [4];
GASSTV_CondatVu.order_filt = [5];

name_method = "GASSTV_CondatVu";

get_params_savetext = @(params) ...
    sprintf("l%.2g_%.2g_sig%s_%g_ns%d_fil%d_k%g_r%.2f_stop1e-%d", ...
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

    HSI_noisy = single(HSI_noisy);


    % Parameter axes you want to analyze
    x_name = "lambda2";          % 横軸にしたい params field 名
    y_name = "sigma_l";          % 縦軸にしたい params field 名
    x_list = GASSTV_CondatVu.(x_name)(:);
    y_list = GASSTV_CondatVu.(y_name)(:);



    % Result grids (rows=sigma_sp, cols=lambda_rho_sp)
    vals_MPSNR = nan(numel(x_list), numel(y_list));
    vals_MSSIM = nan(numel(x_list), numel(y_list));

    [~, ~, ~, info_l] = ...
        Create_SpectralGraphLaplacian(HSI_noisy, 4, "med", 10, 5);
    val_sigma_l_med = info_l.val_sigma_l;

    [~, ~, ~, info_l] = ...
        Create_SpectralGraphLaplacian(HSI_noisy, 4, "90", 10, 5);
    val_sigma_l_90 = info_l.val_sigma_l;


    %% Loop over only (sigma_sp, lambda_rho_sp)
    % Keep others fixed as in your experiment
    for i = 1:numel(x_list)
        for j = 1:numel(y_list)

            params = struct();
            params.lambda_rho_sp = GASSTV_CondatVu.lambda_rho_sp(1);
            params.lambda2 = x_list(j);
            params.sigma_sp = GASSTV_CondatVu.sigma_sp{idx_image};
            params.sigma_l = y_list(i);
            params.num_segments = GASSTV_CondatVu.num_segments(1);
            params.order_filt = GASSTV_CondatVu.order_filt(1);
            params.k_lap = GASSTV_CondatVu.k_lap(1);
            params.maxiter = maxiter;
            params.stopcri = stopcri;
            params.rho_radius = rhos{1};

           
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
    figure;
    t = tiledlayout(1,2, "TileSpacing","compact", "Padding","compact");
    
    % Axes values
    x = log10(y_list(:)).';          % omega_sp
    y = log10(x_list(:));      % log10(sigma_sp)
    
    % Reference sigma values (log-scale)
    y_med = log10(val_sigma_l_med);
    y_90  = log10(val_sigma_l_90);
    
    % Font sizes
    fs_axis  = 13;
    fs_label = 14;
    fs_text  = 14;
    
    % ------------------ MPSNR heatmap ------------------
    nexttile;
    imagesc(x, y, vals_MPSNR);
    set(gca, ...
        "YDir", "normal", ...
        "FontSize", fs_axis, ...
        "XColor", "k", ...
        "YColor", "k");
    
    xlabel("M", "FontSize", fs_label, "Color", "k");
    ylabel("log_{10}(\sigma_{\lambda})", "FontSize", fs_label, "Color", "k");
    
    cb1 = colorbar;
    set(cb1, ...
        "Color", "k", ...
        "FontSize", fs_axis, ...
        "Box", "on");
    
    hold on;
    % yline(y_med, "--k", ...
    %     sprintf(" med : log_{10}(%.3g)", val_sigma_l_med), ...
    %     "LineWidth", 3.0, ...
    %     "FontSize", fs_text, ...
    %     "FontWeight", "bold", ...
    %     "LabelHorizontalAlignment", "left", ...
    %     "LabelVerticalAlignment", "bottom");

    yline(y_med, "--k", ...
        " med", ...
        "LineWidth", 3.0, ...
        "FontSize", fs_text, ...
        "FontWeight", "bold", ...
        "LabelHorizontalAlignment", "left", ...
        "LabelVerticalAlignment", "bottom");
    
    yline(y_90, "-.k", ...
        " 90%", ...
        "LineWidth", 3.0, ...
        "FontSize", fs_text, ...
        "FontWeight", "bold", ...
        "LabelHorizontalAlignment", "left", ...
        "LabelVerticalAlignment", "top");
    hold off;
    
    % ------------------ MSSIM heatmap ------------------
    nexttile;
    imagesc(x, y, vals_MSSIM);
    set(gca, ...
        "YDir", "normal", ...
        "FontSize", fs_axis, ...
        "XColor", "k", ...
        "YColor", "k");
    
    xlabel("\omega", "FontSize", fs_label, "Color", "k");
    ylabel("log_{10}(\sigma_{\lambda})", "FontSize", fs_label, "Color", "k");
    
    cb2 = colorbar;
    set(cb2, ...
        "Color", "k", ...
        "FontSize", fs_axis, ...
        "Box", "on");
    
    hold on;
    yline(y_med, "--k", ...
        sprintf(" med : log_{10}(%.3g)", val_sigma_l_med), ...
        "LineWidth", 3.0, ...
        "FontSize", fs_text, ...
        "FontWeight", "bold", ...
        "LabelHorizontalAlignment", "left", ...
        "LabelVerticalAlignment", "bottom");
    
    yline(y_90, "-.k", ...
        sprintf(" 90%% : log_{10}(%.3g)", val_sigma_l_90), ...
        "LineWidth", 3.0, ...
        "FontSize", fs_text, ...
        "FontWeight", "bold", ...
        "LabelHorizontalAlignment", "left", ...
        "LabelVerticalAlignment", "top");
    hold off;


    %% ===== Output figure (EPS + PNG) =====
    if exist("is_output", "var")

        save_figure_folder = fullfile(...
            dir_save_comp_folder, ...
            "GASSTV_Defence", ...
            "ParamAnal");
    
        % Ensure save directory exists
        if ~exist(save_figure_folder, "dir")
            mkdir(save_figure_folder);
        end
    
        % File base name (paper-friendly)
        fig_base = sprintf( ...
            "%s_GASSTV_g%g_ps%g_pt%g_pd%g_heatmap", ...
            image, ...
            deg.gaussian_sigma, deg.sparse_rate, ...
            deg.stripe_rate, deg.deadline_rate ...
        );
    
        % Full paths (string)
        fig_path_png = fullfile(save_figure_folder, fig_base + ".png");
    
        % ---------- Common figure settings ----------
        set(gcf, "Color", "w");
    
        % ---------- PNG (raster, for preview / slides) ----------
        % Use higher DPI so it looks good in slides
        print(gcf, fig_path_png, "-dpng", "-r300");
    
        fprintf("[Saved PNG] %s\n", fig_path_png);
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
