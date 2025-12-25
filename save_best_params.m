clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

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
% idc_noise_conditions = [7, 6];
idc_noise_conditions = 7;
% idc_noise_conditions = 1:7;
% idc_noise_conditions = 1;
% idc_noise_conditions = [5, 2];

images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

% idc_images = 1:numel(images);
idc_images = 1:2;
% idc_images = 3;
% idc_images = 1;


%% Setting common parameters
% rhos = {0.93, 0.95, 0.98};
rhos = {0.95};
% rhos = {0.98};

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
% GASSTV_CondatVu
% params_GASSTV.sigma_sp = [0.01, 0.1, 1];
GASSTV_CondatVu.sigma_sp = {"90", "med"};
GASSTV_CondatVu.sigma_l = {"med"};
GASSTV_CondatVu.lambda_rho_sp = [0.9, 1, 1.1];
GASSTV_CondatVu.lambda2 = [0.1, 1, 10];
% GASSTV_CondatVu.num_segments = [3, 5, 10];
GASSTV_CondatVu.num_segments = [4];
GASSTV_CondatVu.order_filt = [5, 1];
% GASSTV_CondatVu.order_filt = [5];
methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", false ...
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

% GASSTV_CondatVu_OG
% params_GASSTV.sigma_sp = [0.01, 0.1, 1];
params_GASSTV.sigma_sp = {"90", "med"};
params_GASSTV.sigma_l = {"90", "med"};
params_GASSTV.lambda_rho_sp = [0.9, 1, 1.1];
% params_GASSTV.sigma_s = [0.01, 1];
% params_GASSTV.lambda1 = [0.01, 0.1, 1];
% params_GASSTV.lambda1 = [0.01];
params_GASSTV.lambda2 = [0.1, 1, 10];
params_GASSTV.num_segments = [3, 5, 10];
% params_GASSTV.num_segments = [5];
% params_GASSTV.order_filt = [1, 3];
params_GASSTV.order_filt = [5];
methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu_OG_Med", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{params_GASSTV.lambda_rho_sp, params_GASSTV.lambda2, params_GASSTV.sigma_sp, params_GASSTV.sigma_l, ...
        params_GASSTV.num_segments, params_GASSTV.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% SSTV
methods_info(end+1) = struct( ...
    "name", "SSTV", ...
    "param_names", {{"maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("r%.2f_stop1e-%d", params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% HSSTV_L1
% params_HSSTV.omega = {0.05};
params_HSSTV.omega = {0.01, 0.03, 0.05};

methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "param_names", {{"L", "omega", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{{"L1"}, params_HSSTV.omega, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "param_names", {{"L", "omega", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{{"L12"}, params_HSSTV.omega, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% l0l1HTV
params_l0l1HTV.stepsize_reduction = 0.999;
% params_l0l1HTV.stepsize_reduction = 0.9999;
% params_l0l1HTV.stepsize_reduction = 1;

params_l0l1HTV.L10ball_th = {0.02, 0.03};


methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "param_names", {{"L10ball_th", "stepsize_reduction", ...
        "maxiter", "stopcri", "rho_radius"}}, ...
"params", {{params_l0l1HTV.L10ball_th, params_l0l1HTV.stepsize_reduction, ...
        maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("sr%.5g_th%.2f_r%.2f_maxiter%d", ...
            params.stepsize_reduction, params.L10ball_th, params.rho_radius, maxiter), ...
    "enable", false ...
);

% LRTDTV
params_LRTDTV.tau = 1;
% params_LRTDTV.lambda = 100;
% params_LRTDTV.lambda_param = {10, 15, 20, 25, sqrt(hsi.n1*hsi.n2)};
% params_LRTDTV.lambda = cellfun(@(x) 100 * x / sqrt(hsi.n1*hsi.n2), LRTDTV_lambda_param);
params_LRTDTV.lambda_param = {10, 15, 20, 25, sqrt(100*100)};
params_LRTDTV.lambda = cellfun(@(x) 100 * x / sqrt(100*100), params_LRTDTV.lambda_param);
% params_LRTDTV.rank = {[51,51,10]};
% params_LRTDTV.rank = {[hsi.n1*0.8, hsi.n2*0.8, 10], [51, 51, 10]};
params_LRTDTV.rank = {[100*0.8, 100*0.8, 10], [51, 51, 10]};

methods_info(end+1) = struct( ...
    "name", "LRTDTV", ...
    "param_names", {{"tau", "lambda", "rank"}}, ...
    "params", {{params_LRTDTV.tau, params_LRTDTV.lambda, params_LRTDTV.rank}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.4g_r%d_stop1e-%d", params.lambda, params.rank(1), stopcri_idx), ...
    "enable", false ...
);

% TPTV
params_TPTV.Rank = {[7,7,5]};
params_TPTV.initial_rank = {2};
params_TPTV.maxIter = {50, 100};
% params_TPTV.lambda = 4e-3*sqrt(hsi.n1*hsi.n2);
params_TPTV.lambdas = {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2};

methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "param_names", {{"Rank", "initial_rank", "maxIter", "lambda"}}, ...
    "params", {{params_TPTV.Rank, params_TPTV.initial_rank, params_TPTV.maxIter, params_TPTV.lambdas}}, ...
    "get_params_savetext", @(params) ...
        sprintf("maxiter%d_l%.4g", params.maxIter, params.lambda), ...
    "enable", false ...
);

% FastHyMix
% params_FastHyMix.k_subspace = {4, 6, 8, 10, 12};
params_FastHyMix.k_subspace = {4, 8, 12};

methods_info(end+1) = struct( ...
    "name", "FastHyMix", ...
    "param_names", {{"k_subspace"}}, ...
    "params", {{params_FastHyMix.k_subspace}}, ...
    "get_params_savetext", @(params) ...
        sprintf("sub%d", params.k_subspace), ...
    "enable", false ...
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

summary_metrics_table = struct2table(summary_metrics);

save(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
save(fullfile(dir_result_folder, "summary_metrics.mat"), "summary_metrics", "summary_metrics_table");
writetable(summary_metrics_table, fullfile(dir_result_folder, "summary_metrics.csv"));

fprintf("best param: %s\n", best_params_savetext);
fprintf("MPSNR: %#.4g\n", val_mpsnr_max);


end
end

end