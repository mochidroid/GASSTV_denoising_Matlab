clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Switching operator
is_plot_metrics = 1;
is_show_HSI = 1;
diff_magnification = 7;

is_plot_psnr_and_ssim_per_band = 1;
% is_output_csv = 1;


%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.1,   0,     0,     0,    0.01},  ... % g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % g0.1 ps0.05 pt0.05 pd0.01
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = 7;

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

% images = {...
%     "JasperRidge64", ...
%     "PaviaU64", ...
% };

idc_images = 1:numel(images);
% idc_images = 2;


%% Setting common parameters
rho_radius = 0.95;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;
% maxiter = 5;

load("dir_save_comp_folder.mat", "dir_save_comp_folder");


%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L1
methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L12 
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "line_style", "--", ...
    "enable", true ...
);

% l0l1HTV
methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "line_style", "--", ...
    "enable", true ...
);

% LRTDTV
methods_info(end+1) = struct( ...
    "name", "LRTDTV", ...
    "line_style", "--", ...
    "enable", false ...
);

% TPTV
methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "line_style", "--", ...
    "enable", true ...
);

% FastHyMix
methods_info(end+1) = struct( ...
    "name", "FastHyMix", ...
    "line_style", "--", ...
    "enable", false ...
);

% GASSTV_OraGuide
methods_info(end+1) = struct( ...
    "name", "GASSTV_GLR_Const_OraGuide", ...
    "line_style", "-", ...
    "enable", false ...
);


% GASSTV_NL_OraG
methods_info(end+1) = struct( ...
    "name", "GASSTV_NL_OraG", ...
    "line_style", "-", ...
    "enable", true ...
);


methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Loading results
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Selecting noise condition
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};


[HSI_clean, hsi] = Load_HSI(image);
noise_seed = "default";
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


% Calculation noisy results
val_mpsnr_noisy  = calc_MPSNR(HSI_noisy, HSI_clean);
val_mssim_noisy  = calc_MSSIM(HSI_noisy, HSI_clean);
val_sam_noisy    = calc_SAM(HSI_noisy, HSI_clean);


fprintf("~~~ SETTINGS ~~~\n");
fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n", deg.deadline_rate);


%% Loading results
for idx_method = 1:num_methods
    name_method = methods_info(idx_method).name;

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
        name_method ...
    );

    load(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
    methods_info(idx_method).get_params_savetext = best_params_savetext;

        load(fullfile(dir_result_folder, append(best_params_savetext, ".mat")), ...
            "HSI_restored", "val_mpsnr", "val_mssim", "val_sam");

    methods_info(idx_method).HSI_restored = HSI_restored;
    methods_info(idx_method).val_mpsnr = val_mpsnr;
    methods_info(idx_method).val_mssim = val_mssim;
    methods_info(idx_method).val_sam = val_sam;
end


%% Plotting metrics
if exist("is_plot_metrics", "var") && is_plot_metrics == 1
    % Initialization
    name_length_max = max(strlength([methods_info.name]));
    name_params_length_max = max(strlength([methods_info.get_params_savetext]));

    fprintf("~~~ RESULTS ~~~\n");
    fprintf("%s  \t MPSNR\t MSSIM\t SAM\n", blanks(name_length_max + name_params_length_max + 2));

    % Plotting noisy results
    fprintf("%s: \t %#.4g\t %#.4g\t %#.4g\n", ...
        append("noisy", blanks(name_length_max + name_params_length_max - 3)), ...
        val_mpsnr_noisy, val_mssim_noisy, val_sam_noisy);

    % Plotting restored results
    for idx_method = 1:num_methods
        name_method         = methods_info(idx_method).name;
        params_save_text    = methods_info(idx_method).get_params_savetext;
        val_mpsnr           = methods_info(idx_method).val_mpsnr;
        val_mssim           = methods_info(idx_method).val_mssim;
        val_sam             = methods_info(idx_method).val_sam;

        fprintf("%s(%s): \t %#.4g\t %#.4g\t %#.4g\n", ...
            append(name_method, blanks(name_length_max - strlength(name_method))), ...
            append(params_save_text, blanks(name_params_length_max - strlength(params_save_text))), ...
            val_mpsnr, val_mssim, val_sam);
    end
end


%% Showing HSI
if exist("is_show_HSI", "var") && is_show_HSI == 1
    cat_HSI = cat(2, HSI_clean, HSI_noisy);

    for idx_method = 1:num_methods
        cat_HSI = cat(2, cat_HSI, methods_info(idx_method).HSI_restored);
    end

    cat_diff = abs(repmat(HSI_clean, [1, num_methods+2, 1]) - cat_HSI) * diff_magnification;

    implay(cat(1, cat_HSI, cat_diff));
end


%% Plotting psnr and ssim per band
if exist("is_plot_psnr_and_ssim_per_band", "var") && is_plot_psnr_and_ssim_per_band == 1
    % Calculating psnr and ssim per band
    vals_psnr_per_band = zeros(num_methods, hsi.n3);
    vals_ssim_per_band = zeros(num_methods, hsi.n3);

    for idx_method = 1:num_methods
        [psnr_per_band, ssim_per_band] = ...
            calc_PSNR_SSIM_per_band(methods_info(idx_method).HSI_restored, HSI_clean);
        vals_psnr_per_band(idx_method, :) = psnr_per_band;
        vals_ssim_per_band(idx_method, :) = ssim_per_band;
    end

    % Plotting results
    label_style = {"FontSize", 15};
    
    figure()
    hold on
    for idx_method = 1:num_methods
        name_method = methods_info(idx_method).name;
        lineStyle_method = methods_info(idx_method).line_style;
        plot_style = {"LineWidth", 2.0, "LineStyle", lineStyle_method, "DisplayName", name_method};
        plot(vals_psnr_per_band(idx_method,:), plot_style{:})
    end
    hold off
    title("psnr per band")
    xlabel("band", label_style{:})
    ylabel("psnr", label_style{:})
    
    figure()
    hold on
    for idx_method = 1:num_methods
        name_method = methods_info(idx_method).name;
        lineStyle_method = methods_info(idx_method).line_style;
        plot_style = {"LineWidth", 2.0, "LineStyle", lineStyle_method, "DisplayName", name_method};
        plot(vals_ssim_per_band(idx_method,:), plot_style{:})
    end
    hold off
    title("ssim per band")
    xlabel("band", label_style{:})
    ylabel("ssim", label_style{:})

end


%% Outputting csv
if exist("is_output_csv", "var") && ~isempty(is_output_csv)
    %initial
    dir_condition_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)) ...
    );
    fmt4s = @(x) round(x,4,"significant");

    % noisy
    tbl_noisy = table(...
        "noisy", "", fmt4s(val_mpsnr_noisy), fmt4s(val_mssim_noisy), fmt4s(val_sam_noisy), ...
        'VariableNames', ["method","params","mpsnr","mssim","sam"] ...
    );
    tbl_list = tbl_noisy;

    %each method
    for idx_method = 1:num_methods
        name_method = methods_info(idx_method).name;
        name_params = methods_info(idx_method).get_params_savetext;
        val_mpsnr   = fmt4s(methods_info(idx_method).val_mpsnr);
        val_mssim   = fmt4s(methods_info(idx_method).val_mssim);
        val_sam     = fmt4s(methods_info(idx_method).val_sam);
        tbl_method = table(...
            name_method, name_params, val_mpsnr, val_mssim, val_sam, ...
            'VariableNames', ["method","params","mpsnr","mssim","sam"] ...
        );
        tbl_list = [tbl_list; tbl_method];
    end
    writetable(tbl_list, fullfile(dir_condition_folder, 'summary_metrics.csv'));
end

end
end
