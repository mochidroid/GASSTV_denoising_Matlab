clear;
close all;

addpath(genpath('sub_functions'))


%% Selecting conditions
noise_conditions = { ...
    {0.1,   0,      0,      0},     ... % g0.1 ps0 pt0
    {0,     0,      0.05,   0.5},   ... % g0 ps0 pt0.05
    {0.05,  0.05,   0,      0},     ... % g0.05 ps0.05 pt0
    {0.1,   0.05,   0,      0},     ... % g0.1 ps0.05 pt0
    {0.05,  0,      0.05,   0.5},   ... % g0.05 ps0 pt0.05
    {0.1,   0,      0.05,   0.5},   ... % g0.1 ps0 pt0.05
    {0.05,  0.05,   0.05,   0.5},   ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,   0.05,   0.5},   ... % g0.1 ps0.05 pt0.05
};

idc_noise_conditions = 1:2;


images = {...
    'JasperRidge', ...
    'PaviaU120', ...
    'Beltsville', ...
};

% idc_images = 1:numel(images);
idc_images = 1:2;


for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Selecting noise condition
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
image = images{idx_image};

[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


diff_magnification = 7;


%% Selecting methods
name_method = 'TPTV';

fprintf('\n~~~ SETTINGS ~~~\n');
fprintf('Method: %s\n', name_method);
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);
fprintf('Gaussian sigma: %g\n', deg.gaussian_sigma);
fprintf('Sparse rate: %g\n', deg.sparse_rate);
fprintf('Stripe rate: %g\n', deg.stripe_rate);
fprintf('Stripe intensity: %g\n', deg.stripe_intensity);


%% Setting parameters
rho = 0.95;
stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;
maxiter = 20000;

% TPTV
TPTV_Rank = {[7,7,5]};
TPTV_initial_rank = {2};
TPTV_maxIter = {50, 100};
% TPTV_lambda = 4e-3*sqrt(hsi.n1*hsi.n2);
TPTV_lambdas = {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2};

params_tmp = {TPTV_Rank, TPTV_initial_rank, TPTV_maxIter, TPTV_lambdas}; % TPTV
name_params = {'Rank', 'initial_rank', 'maxIter', 'lambda'}; % TPTV

[params_comb, num_params_comb] = ParamsList2Comb(params_tmp);


%% Computing savetext for loading
names_params_savetext = strings(num_params_comb, 1);

for idx_params_comb = 1:num_params_comb
    params = struct();
    for idx_params = 1:numel(name_params)
        % Assigning parameters to the structure
        params.(name_params{idx_params}) = params_comb{idx_params_comb}{idx_params};
    end
    
    names_params_savetext(idx_params_comb) = ...
        append('maxiter', num2str(params.maxIter), '_l', num2str(params.lambda)); % TPTV
end

names_params_savetext_max = max(strlength(names_params_savetext), [], 'all');

fprintf('~~~ RESULTS ~~~\n');
fprintf('   %s   \t MPSNR\t MSSIM\t SAM\n', blanks(names_params_savetext_max));


%% Initialiging best param
name_params_savetext_best_mpsnr   = struct();

best_mpsnr  = 0;


%% Loading restored HS images
for idx_params_comb = 1:num_params_comb
    load_folder_name = append(...
        './result/' , ...
        'denoising_', image, '/', ...
        'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
            '_pt', num2str(deg.stripe_rate), '/', ...
        name_method, '/', ...
        names_params_savetext(idx_params_comb), '/' ...   
    );
    
    load(append(load_folder_name, 'metric_vals.mat'), ...
        'mpsnr', 'mssim', 'sam' ...
    );
    
    fprintf('%2d. %s: \t %#.4g\t %#.4g\t %#.4g\n', ...
        idx_params_comb, append(names_params_savetext(idx_params_comb), ...
        blanks(names_params_savetext_max - strlength(names_params_savetext(idx_params_comb)))), ...
        mpsnr, mssim, sam);


    if best_mpsnr < mpsnr
        best_mpsnr = mpsnr;
    
        name_params_savetext_best_mpsnr = names_params_savetext(idx_params_comb);
    end
end


%% Saving best mpsnr restored HS images
load_best_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method, '/', ...
    name_params_savetext_best_mpsnr, '/' ...   
);

load(append(load_best_folder_name, 'image_result.mat'));
load(append(load_best_folder_name, 'metric_vals.mat'));
load(append(load_best_folder_name, 'other_result.mat'));


fprintf('~~~ BEST MPSNR RESULTS ~~~\n');
fprintf('Parameter: %s\n', name_params_savetext_best_mpsnr);
fprintf('MPSNR: %#.4g, MSSIM: %#.4g, SAM: %#.4g\n', mpsnr, mssim, sam);


is_rmdir_best_params = 1;
if exist('is_redir_best_params', 'var') && is_rmdir_best_params
    save_best_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method, '/', ...
    '/best_params/' ...   
    );

    if exist(save_best_folder_name, 'dir')
        rmdir(save_best_folder_name, 's');
    end
end

save_folder_name_for_best = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method, '/', ...
    '/best_params/', ...
    name_params_savetext_best_mpsnr, '/' ...   
);

mkdir(save_folder_name_for_best);

save(append(save_folder_name_for_best, 'image_result.mat'), ...
    'HSI_clean', 'HSI_noisy', 'hsi', 'deg', 'image', ...
    'HSI_restored', 'removed_noise', ...
    '-v7.3', '-nocompression' ...
);

save(append(save_folder_name_for_best, 'metric_vals.mat'), ...
    'mpsnr', 'mssim', 'ergas', 'sam', ...
    'psnr_per_band', 'ssim_per_band', ...
    '-v7.3', '-nocompression' ...
);

save(append(save_folder_name_for_best, 'other_result.mat'), ...
    'other_result', ...
    '-v7.3', '-nocompression' ...
);


end
end
