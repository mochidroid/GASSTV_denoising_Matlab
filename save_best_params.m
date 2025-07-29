clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

fprintf("******* initium *******\n");

%% Selecting conditions
noise_conditions = { ...
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.1,   0.1,   0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0,     0,     0.05,  0.5,  0},     ... % g0 ps0 pt0.05
    {0.05,  0.05,  0,     0,    0},     ... % g0.05 ps0.05 pt0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.3,  0},     ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,  0.05,  0.3,  0},     ... % g0.1 ps0.05 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05
    {0.05,  0.05,  0.05,  0.3,  0.001}, ... % g0.05 ps0.05 pt0.05 pd0.001
    {0.1,   0.05,  0.05,  0.3,  0.001}, ... % g0.1 ps0.05 pt0.05 pd0.001
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
% idc_noise_conditions = [5:6, 11:14];
% idc_noise_conditions = 5:6;
idc_noise_conditions = 11:12;

images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

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


%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "param_names", {{"maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("r%.2f_stop1e-%d", params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% HSSTV_L1
% HSSTV_omega = {0.05};
HSSTV_omega = {0.01, 0.03, 0.05};

methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "param_names", {{"L", "omega", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{{"L1"}, HSSTV_omega, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "param_names", {{"L", "omega", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{{"L12"}, HSSTV_omega, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% l0l1HTV
l0l1HTV_stepsize_reduction = 0.999;
% l0l1HTV_stepsize_reduction = 0.9999;
% l0l1HTV_stepsize_reduction = 1;

l0l1HTV_L10ball_th = {0.02, 0.03};


methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "param_names", {{"L10ball_th", "stepsize_reduction", ...
        "maxiter", "stopcri", "rho_radius"}}, ...
"params", {{l0l1HTV_L10ball_th, l0l1HTV_stepsize_reduction, ...
        maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("sr%.5g_th%.2f_r%.2f_maxiter%d", ...
            params.stepsize_reduction, params.L10ball_th, params.rho_radius, maxiter), ...
    "enable", false ...
);

% LRTDTV
LRTDTV_tau = 1;
% LRTDTV_lambda = 100;
% LRTDTV_lambda_param = {10, 15, 20, 25, sqrt(hsi.n1*hsi.n2)};
% LRTDTV_lambda = cellfun(@(x) 100 * x / sqrt(hsi.n1*hsi.n2), LRTDTV_lambda_param);
LRTDTV_lambda_param = {10, 15, 20, 25, sqrt(100*100)};
LRTDTV_lambda = cellfun(@(x) 100 * x / sqrt(100*100), LRTDTV_lambda_param);
% LRTDTV_rank = {[51,51,10]};
% LRTDTV_rank = {[hsi.n1*0.8, hsi.n2*0.8, 10], [51, 51, 10]};
LRTDTV_rank = {[100*0.8, 100*0.8, 10], [51, 51, 10]};

methods_info(end+1) = struct( ...
    "name", "LRTDTV", ...
    "param_names", {{"tau", "lambda", "rank"}}, ...
    "params", {{LRTDTV_tau, LRTDTV_lambda, LRTDTV_rank}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.4g_r%d_stop1e-%d", params.lambda, params.rank(1), stopcri_idx), ...
    "enable", false ...
);

% TPTV
TPTV_Rank = {[7,7,5]};
TPTV_initial_rank = {2};
TPTV_maxIter = {50, 100};
% TPTV_lambda = 4e-3*sqrt(hsi.n1*hsi.n2);
TPTV_lambdas = {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2};

methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "param_names", {{"Rank", "initial_rank", "maxIter", "lambda"}}, ...
    "params", {{TPTV_Rank, TPTV_initial_rank, TPTV_maxIter, TPTV_lambdas}}, ...
    "get_params_savetext", @(params) ...
        sprintf("maxiter%d_l%.4g", params.maxIter, params.lambda), ...
    "enable", false ...
);

% GASSTV
GASSTV_sigma_sp = [0.01, 0.1, 1];
GASSTV_sigma_s = [0.01, 1];
GASSTV_lambda1 = [0.01, 0.1, 1];
% GASSTV_lambda1 = [0.01];
GASSTV_lambda2 = [0.01, 0.1, 1];
% GASSTV_num_segments = [3, 5, 10];
GASSTV_num_segments = [3];

methods_info(end+1) = struct( ...
    "name", "GASSTV", ...
    "param_names", {{"L", "lambda1", "lambda2", "sigma_sp", "sigma_s", ...
        "num_segments", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{{"L1"}, GASSTV_lambda1, GASSTV_lambda2, GASSTV_sigma_sp, GASSTV_sigma_s, ...
        GASSTV_num_segments, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sp%.2g_s%.2g_ns%d_r%.2f_stop1e-%d", ...
        params.lambda1, params.lambda2, params.sigma_sp, params.sigma_s, ...
        params.num_segments, params.rho_radius, stopcri_idx), ...
    "enable", false ...
);


% GASSTV_Const
GASSTV_sigma_sp = [0.01, 0.05, 0.1];
GASSTV_sigma_l = [0.01, 0.05, 0.1];
GASSTV_lambda_rho_sp = [0.8, 0.9, 1, 1.1];
% GASSTV_lambda1 = [0.01];
GASSTV_lambda_rho_l = [0.8, 0.9, 1, 1.1];
% GASSTV_num_segments = [3, 5, 10];
GASSTV_num_segments = [3, 5];

methods_info(end+1) = struct( ...
    "name", "GASSTV_Const", ...
    "param_names", {{"sigma_sp", "sigma_l", "lambda_rho_sp", "lambda_rho_l", ...
        "num_segments", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_sigma_sp, GASSTV_sigma_l, GASSTV_lambda_rho_sp, GASSTV_lambda_rho_l, ...
        GASSTV_num_segments, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("lr%.2g_%.2g_sp%.2g_s%.2g_ns%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda_rho_l, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.rho_radius, stopcri_idx), ...
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

for idx_params_comb = 1:num_params_comb
    params = struct();
    for idx_params = 1:numel(params_name)
        % Assigning parameters to the structure
        params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
    end

    names_params_savetext(idx_params_comb) = ...
        methods_info(idx_method).get_params_savetext(params);

end

names_params_savetext_max = max(strlength(names_params_savetext), [], "all");


fprintf("~~~ RESULTS ~~~\n");
fprintf("   %s   \t MPSNR\t MSSIM\t SAM\n", blanks(names_params_savetext_max));


%% Calculating best mpsnr
% Initialization
best_params_savetext = strings(1);
val_mpsnr_max  = 0;

% Calculating best mpsnr
for idx_params_comb = 1:num_params_comb
    name_params_savetext = names_params_savetext(idx_params_comb);

    load(fullfile(dir_result_folder, append(name_params_savetext, ".mat")), ...
        "val_mpsnr", "val_mssim", "val_sam");
    
    fprintf("%2d. %s: \t %#.4g\t %#.4g\t %#.4g\n", ...
        idx_params_comb, append(name_params_savetext, ...
        blanks(names_params_savetext_max - strlength(name_params_savetext))), ...
        val_mpsnr, val_mssim, val_sam);


    if val_mpsnr_max < val_mpsnr
        val_mpsnr_max = val_mpsnr;
        best_params_savetext = name_params_savetext;
    end

end


save(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");

fprintf("best param: %s\n", best_params_savetext);
fprintf("MPSNR: %#.4g\n", val_mpsnr_max);


end
end

end