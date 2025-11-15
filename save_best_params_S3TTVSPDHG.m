clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

fprintf("******* initium *******\n");

%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.1 pt0.05 pd0
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = 2;

images = {...
    "JasperRidge64", ...
    "PaviaU64", ...
};

% images = {...
%     "JasperRidge", ...
%     "PaviaU", ...
%     "Beltsville", ...
% };

% idc_images = 1:numel(images);
idc_images = 1;


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

fmt4s = @(x) round(x,4,"significant");


%% Setting each methods info
% S3TTV_SPDHG
params_S3TTV.blocksize = {[8,8]};
% params_S3TTV.blocksize = {[10,10]};
% params_S3TTV.prob_patch = {1};
params_S3TTV.prob_patch = {0.25, 0.1, 0.5};
% params_S3TTV.prob_patch = {0.25, 0.1, 0.5, 1};
params_S3TTV.eta = {0};
% params_S3TTV.eta = {0.9, 0};
params_S3TTV.clip_c = {1};
% params_S3TTV.clip_c = {2};

% methods_info(1) = struct( ...
%     "name", "S3TTV_SPDHG_v5", ...
%     "param_names", {{"blocksize", "eta", "prob_patch", "maxiter", "stopcri", "rho_radius"}}, ...
%     "params", {{params_S3TTV.blocksize, params_S3TTV.eta, params_S3TTV.prob_patch, maxiter, stopcri, rhos}}, ...
%     "get_params_savetext", @(params) ...
%         sprintf("bl%d_p%g_e%g_r%.2f_stop1e-%d", ...
%         params.blocksize(1), params.prob_patch, params.eta, params.rho_radius, stopcri_idx), ...
%     "enable", true ...
% );

methods_info(1) = struct( ...
    "name", "S3TTV_SPDHG_v5", ...
    "param_names", {{"blocksize", "eta", "clip_c", "prob_patch", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{params_S3TTV.blocksize, params_S3TTV.eta, params_S3TTV.clip_c, params_S3TTV.prob_patch, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("bl%d_p%g_e%g_c%g_r%.2f_stop1e-%d", ...
        params.blocksize(1), params.prob_patch, params.eta, params.clip_c, params.rho_radius, stopcri_idx), ...
    "enable", true ...
);


methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);

load("dir_save_comp_folder.mat", "dir_save_comp_folder");


%% Loading each method
for idx_method = 1:num_methods
name_method = methods_info(idx_method).name;
params_name = methods_info(idx_method).param_names;
params_cell = methods_info(idx_method).params;

[params_comb, num_params_comb] = ParamsList2Comb(params_cell);


%% Loading each condition
for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Generating observation
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};

% [HSI_clean, hsi] = Load_HSI(image);
% noise_seed = "default";
% [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);
% 
% HSI_clean = single(HSI_clean);
% HSI_noisy = single(HSI_noisy);


dir_result_folder = fullfile(...
    dir_save_comp_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);


fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Method: %s\n", name_method);
fprintf("Image: %s \n", image);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n", deg.deadline_rate);


%% Computing savetext for loading
names_params_savetext = strings(num_params_comb, 1);
clear summary_metrics

for idx_params_comb = 1:num_params_comb
    params = struct();
    for idx_params = 1:numel(params_name)
        % Assigning parameters to the structure
        params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
    end

    if idx_params_comb == 1
        summary_metrics(1) = params;
    else
        summary_metrics(end+1) = params;
    end

    names_params_savetext(idx_params_comb) = ...
        methods_info(idx_method).get_params_savetext(params);

end

names_params_savetext_max = max(strlength(names_params_savetext), [], "all");


fprintf("~~~ RESULTS ~~~\n");
fprintf("   %s   \t MPSNR\t MSSIM\t SAM\n", blanks(names_params_savetext_max));


%% Calculating best mpsnr
% Calculating best mpsnr
for idx_params_comb = 1:num_params_comb
    name_params_savetext = names_params_savetext(idx_params_comb);

    load(fullfile(dir_result_folder, append(name_params_savetext, ".mat")), ...
        "val_mpsnr", "val_mssim", "val_sam");
    
    fprintf("%2d. %s: \t %#.4g\t %#.4g\t %#.4g\n", ...
        idx_params_comb, append(name_params_savetext, ...
        blanks(names_params_savetext_max - strlength(name_params_savetext))), ...
        val_mpsnr, val_mssim, val_sam);

    summary_metrics(idx_params_comb).mpsnr = fmt4s(val_mpsnr);
    summary_metrics(idx_params_comb).mssim = fmt4s(val_mssim);
    summary_metrics(idx_params_comb).sam   = fmt4s(val_sam);

end


val_mpsnr_vec = [summary_metrics.mpsnr];
[val_mpsnr_max, best_param_index] = max(val_mpsnr_vec);

best_params_savetext = names_params_savetext(best_param_index);

save(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
save(fullfile(dir_result_folder, "summary_metrics.mat"), "summary_metrics");

fprintf("best param: %s\n", best_params_savetext);
fprintf("MPSNR: %#.4g\n", val_mpsnr_max);


end
end

end