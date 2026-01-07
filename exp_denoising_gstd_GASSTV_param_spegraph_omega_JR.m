clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")
addpath("methods")

fprintf("******* initium *******\n");

%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % g0.05 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
    {0.1,   0,     0,     0,    0.01},  ... % g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % g0.1 ps0.05 pt0.05 pd0.01
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
% idc_noise_conditions = [5];
idc_noise_conditions = [7];

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

% images = {...
%     "JasperRidge64", ...
%     "PaviaU64", ...
% };

% idc_images = 1:numel(images);
idc_images = 1;


idx_exp = 0;
total_exp = length(idc_noise_conditions) * length(idc_images);

load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder");


%% Setting common parameters
% rhos = {0.93, 0.95, 0.98};
rhos = {0.95};

% epsilon_rho = 0.01;

% % epsilon = epsilon_rho * sqrt(hsi.N * (1 - deg.sparse_rate));
% epsilon = rhos * deg.gaussian_sigma * sqrt(hsi.N * (1 - deg.sparse_rate));
% alpha = rhos * (0.5 * hsi.N * deg.sparse_rate);
% beta = rhos * hsi.N * (1 - deg.sparse_rate) ...
%     * deg.stripe_rate * deg.stripe_intensity / 2;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;
% maxiter = 5;


%% Setting each methods info
% GASSTV_CondatVu
% GASSTV_CondatVu_OG_Med.sigma_sp = [0.01, 0.1, 1];
% GASSTV_CondatVu.sigma_sp = {"90", "med"};
GASSTV_CondatVu.sigma_sp = {"med"};
% GASSTV_CondatVu.lambda_rho_sp = [0.9, 1, 1.1];
GASSTV_CondatVu.lambda_rho_sp = [0.9];

% GASSTV_CondatVu.sigma_l = {"med"};
GASSTV_CondatVu.sigma_l = {0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1};
% GASSTV_CondatVu.lambda2 = [1];
% GASSTV_CondatVu.lambda2 = [0.01, 0.05, 0.1, 0.5, 1, 5, 10];
% GASSTV_CondatVu.lambda2 = [0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500];
GASSTV_CondatVu.lambda2 = [5, 10, 50, 100];
% GASSTV_CondatVu.k_lap = [10];
GASSTV_CondatVu.k_lap = [1000];
% GASSTV_CondatVu.k_lap = [2, 5, 8, 10, 15, 20, 50, 100];
GASSTV_CondatVu.num_segments = [4];
% GASSTV_CondatVu.num_segments = [2, 3, 4, 5, 6, 7, 8];
% GASSTV_CondatVu.order_filt = [5, 1];
GASSTV_CondatVu.order_filt = [5];
methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
        select_func_GASSTV_for_denoising_CondatVu(HSI_clean, HSI_noisy, params, deg), ...
    "param_names", {{"lambda_rho_sp", "lambda2", ...
        "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "k_lap", ...
        "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, ...
        GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, GASSTV_CondatVu.k_lap, ...
        maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%g_ns%d_fil%d_k%g_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.k_lap, ...
        params.rho_radius, stopcri_idx), ...
    "enable", true ...
);


% GASSTV_CondatVu_OG
% GASSTV_CondatVu_OG.sigma_sp = [0.01, 0.1, 1];
GASSTV_CondatVu_OG.sigma_sp = {"90", "med"};
% GASSTV_CondatVu_OG.sigma_l = {"90", "med"};
GASSTV_CondatVu_OG.sigma_l = {"med"};
GASSTV_CondatVu_OG.lambda_rho_sp = [0.9, 1, 1.1];
% GASSTV_CondatVu_OG.sigma_s = [0.01, 1];
% GASSTV_CondatVu_OG.lambda1 = [0.01, 0.1, 1];
% GASSTV_CondatVu_OG.lambda1 = [0.01];
GASSTV_CondatVu_OG.lambda2 = [0.1, 1, 10];
% GASSTV_CondatVu_OG.num_segments = [3, 5, 10];
GASSTV_CondatVu_OG.num_segments = [4];
% GASSTV_CondatVu_OG.order_filt = [1, 3];
GASSTV_CondatVu_OG.order_filt = [1];
methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu_OG", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
        select_func_GASSTV_OG_for_denoising_CondatVu(HSI_clean, HSI_noisy, params, deg), ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_CondatVu_OG.lambda_rho_sp, GASSTV_CondatVu_OG.lambda2, GASSTV_CondatVu_OG.sigma_sp, GASSTV_CondatVu_OG.sigma_l, ...
        GASSTV_CondatVu_OG.num_segments, GASSTV_CondatVu_OG.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);

% GASSTV_CondatVu_OG_Med
% GASSTV_CondatVu_OG_Med.sigma_sp = [0.01, 0.1, 1];
GASSTV_CondatVu_OG_Med.sigma_sp = {"90", "med"};
GASSTV_CondatVu_OG_Med.sigma_l = {"med"};
GASSTV_CondatVu_OG_Med.lambda_rho_sp = [0.9, 1, 1.1];
% GASSTV_CondatVu_OG_Med.sigma_s = [0.01, 1];
% GASSTV_CondatVu_OG_Med.lambda1 = [0.01, 0.1, 1];
% GASSTV_CondatVu_OG_Med.lambda1 = [0.01];
GASSTV_CondatVu_OG_Med.lambda2 = [0.1, 1, 10];
% GASSTV_CondatVu_OG_Med.num_segments = [3, 5, 10];
GASSTV_CondatVu_OG_Med.num_segments = [4];
% GASSTV_CondatVu_OG_Med.order_filt = [1, 3];
GASSTV_CondatVu_OG_Med.order_filt = [5];
methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu_OG_Med", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
        select_func_GASSTV_OG_for_denoising_CondatVu(HSI_clean, HSI_noisy, params, deg), ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_CondatVu_OG_Med.lambda_rho_sp, GASSTV_CondatVu_OG_Med.lambda2, GASSTV_CondatVu_OG_Med.sigma_sp, GASSTV_CondatVu_OG_Med.sigma_l, ...
        GASSTV_CondatVu_OG_Med.num_segments, GASSTV_CondatVu_OG_Med.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);



methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Running Expt.
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Generating observation
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

idx_exp = idx_exp + 1;

i_method = 0;


%% Running methods
for idx_method = 1:num_methods
name_method = methods_info(idx_method).name;
func_method = methods_info(idx_method).func;
params_name = methods_info(idx_method).param_names;
params_cell = methods_info(idx_method).params;


i_method = i_method + 1;


[params_comb, num_params_comb] = ParamsList2Comb(params_cell);

for idx_params_comb = 1:num_params_comb

params = struct();
for idx_params = 1:numel(params_name)
    % Assigning parameters to the structure
    params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};

    % If rho_radius exists, calculate epsilon, alpha, and beta
    if strcmp(params_name{idx_params}, "rho_radius")
        % Calculating the rate except sparse noise and deadline noise
        rate_except_sd = 1 - deg.sparse_rate - 2*deg.deadline_rate + 2*deg.sparse_rate*deg.deadline_rate;
        params.epsilon = params.rho_radius * deg.gaussian_sigma * sqrt(hsi.N * rate_except_sd);
        params.alpha = params.rho_radius * (0.5 * hsi.N * (deg.sparse_rate - 2*deg.sparse_rate*deg.deadline_rate));
        params.beta = params.rho_radius * hsi.N * rate_except_sd ...
            * deg.stripe_rate * deg.stripe_intensity / 2;
    end
end

name_params_savetext = methods_info(idx_method).get_params_savetext(params);


fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Method: %s\n", name_method);
fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n", deg.deadline_rate);
fprintf("Parameter settings: %s\n", name_params_savetext)
fprintf("Methods: (%d/%d), Cases: (%d/%d), Params:(%d/%d)\n", ...
    i_method, num_methods, idx_exp, total_exp, idx_params_comb, num_params_comb);


dir_save_method_folder = fullfile( ...
    dir_save_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);

dir_save_comp_method_folder = fullfile( ...
    dir_save_comp_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);

file_comp = fullfile(dir_save_comp_method_folder, append(name_params_savetext, ".mat"));

if isfile(file_comp)
    fprintf("[SKIP] Already exists: %s | %s\n", name_method, name_params_savetext);
    continue;  % idx_params_comb ループをスキップ
end


% Run
[HSI_restored, removed_noise, other_result]...
    = func_method(HSI_clean, HSI_noisy, params, deg);


% Plotting results
val_mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
val_mssim  = calc_MSSIM(HSI_restored, HSI_clean);
val_sam    = calc_SAM(HSI_restored, HSI_clean);

fprintf("~~~ RESULTS ~~~\n");
fprintf("MPSNR: %#.4g\n", val_mpsnr);
fprintf("MSSIM: %#.4g\n", val_mssim);
fprintf("SAM  : %#.4g\n", val_sam);

[vals_psnr_per_band, vals_ssim_per_band] = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


% Saving each result
mkdir(dir_save_method_folder);
mkdir(dir_save_comp_method_folder);

save(fullfile(dir_save_method_folder, append(name_params_savetext, ".mat")), ...
    "HSI_clean", "HSI_noisy", "hsi", "deg", "image", ...
    "HSI_restored", "removed_noise", ...
    "val_mpsnr", "val_mssim", "val_sam", ...
    "vals_psnr_per_band", "vals_ssim_per_band", ...
    "params", "other_result", ...
    "-v7.3", "-nocompression" ...
);
save(fullfile(dir_save_comp_method_folder, append(name_params_savetext, ".mat")), ...
    "HSI_restored", "params", "val_mpsnr", "val_mssim", "val_sam", ...
    "-v7.3", "-nocompression" ...
);

close all

end

end
end
end

fprintf("******* finis *******\n");
