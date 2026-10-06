clear
close all;

% addpath('function_S3TTV')
addpath(genpath('sub_functions'))
addpath('func_metrics')
addpath(genpath('methods'))

fprintf('******* initium *******\n');
rng('default')

%% Selecting conditions
images = {...
    'IndianPines', ...
    'Swannee', ...
};

idc_images = 1:numel(images);
% idc_images = 1;


idx_exp = 0;
total_exp = length(idc_images);


for idx_image = idc_images
%% Generating observation
image = images{idx_image};

[HSI_noisy, hsi] = Load_real_HSI(image);
HSI_noisy = single(HSI_noisy);

idx_exp = idx_exp + 1;


%% Selecting methods
names_methods = { ...
    'S3TTV', ...
    'SSTV', ...
    'HSSTV_L1', ...
    'HSSTV_L12', ...
    'l0l1HTV', ...
    'STV', ...
    'SSST', ...
    'LRTDTV', ...
    'TPTV', ...
};


func_methods = {...
    @(params) func_S3TTV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_SSTV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_HSSTV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_HSSTV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_l0l1HTV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_STV_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_SSST_real_for_denoising_by_PsPPDS(HSI_noisy, params), ...
    @(params) func_LRTDTV(HSI_noisy, params), ...
    @(params) func_TPTV_real_for_denoising(HSI_noisy, params), ...
};

idc_methods = 1:numel(names_methods);

i_method = 0;



%% Setting common parameters
switch image
    case 'IndianPines'
        epsilon = 30;
        alpha = 200;
        beta = 100;

    case 'Swannee'
        epsilon = 100;
        alpha = 800;
        beta = 5000;
end


stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;

maxiter = 20000;
% maxiter = 10;

% boundary = 'Neumann';
boundary = {'circulant'};

blocksize = {[10,10]}; % block size of STV-type methods
shiftstep = {[1,1]}; % [1,1] means full overlap



%% Setting unique parameters
% HSSTV
HSSTV_L1_L = {'L1'};
HSSTV_L12_L = {'L12'};
HSSTV_omega = {0.05};


% l0l1HTV
l0l1HTV_stepsize_reduction = 0.999;
% l0l1HTV_stepsize_reduction = 0.9999;
% l0l1HTV_stepsize_reduction = 1;

l0l1HTV_L10ball_th = 0.03;

% L0gradientのパラメータであるthetaを決定する
HSI_mean_3 = mean(HSI_noisy, 3);
diff_mean = (abs(Dv_Neumann(HSI_mean_3)) + abs(Dh_Neumann(HSI_mean_3)))/2;
% diff_mean = (abs(Dv_circulant(HSI_mean_3)) + abs(Dh_circulant(HSI_mean_3)))/2;

diff_mean(diff_mean < l0l1HTV_L10ball_th) = 0; 
% gamma_prime = nnz(diff_mean)/(128*128); % gamma'の導出方法^
% theta = floor(128*128*gamma_prime);
l0l1HTV_eps_L10ball = nnz(diff_mean);


% STV-type methods
STV_tmp_save_name   = append('./tmp_save/', image, '_STV');
SSST_tmp_save_name  = append('./tmp_save/', image, '_SSST');
S3TTV_tmp_save_name = append('./tmp_save/', image, '_S3TTV');


% LRTDTV
LRTDTV_tau = 1;
LRTDTV_lambda = 100;
% LRTDTV_lambda_param = {10, 15, 20, 25, sqrt(hsi.n1*hsi.n2)};
% LRTDTV_lambda = cellfun(@(x) 100 * x / sqrt(hsi.n1*hsi.n2), LRTDTV_lambda_param);
LRTDTV_rank = {[51,51,10]};
% LRTDTV_rank = {[hsi.n1*0.8, hsi.n2*0.8, 10], [51, 51, 10]};


% TPTV
TPTV_Rank = {[7,7,5]};
TPTV_initial_rank = {2};
TPTV_maxIter = {50, 100};
% TPTV_lambda = {4e-3*sqrt(hsi.n1*hsi.n2)};
TPTV_lambdas = {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2, 4e-3*sqrt(hsi.n1*hsi.n2)};


%% Cordinating parameters
params_all = { ...
    {epsilon, alpha, beta, boundary, blocksize, shiftstep, ...
        maxiter, stopcri, S3TTV_tmp_save_name}, ... % S3TTV
    {epsilon, alpha, beta, boundary, maxiter, stopcri}, ... % SSTV
    {epsilon, alpha, beta, boundary, HSSTV_L1_L, HSSTV_omega, maxiter, stopcri}, ... % HSSTV_L1
    {epsilon, alpha, beta, boundary, HSSTV_L12_L, HSSTV_omega, maxiter, stopcri}, ... % HSSTV_L12
    {epsilon, alpha, beta, l0l1HTV_eps_L10ball, l0l1HTV_stepsize_reduction, ...
        boundary, maxiter, stopcri}, ... % l0l1HTV
    {epsilon, alpha, beta, boundary, blocksize, shiftstep, ...
        maxiter, stopcri, STV_tmp_save_name}, ... % STV
    {epsilon, alpha, beta, boundary, blocksize, shiftstep, ...
        maxiter, stopcri, SSST_tmp_save_name}, ... % SSST
    {LRTDTV_tau, LRTDTV_lambda, LRTDTV_rank}, ... % LRTDTV
    {TPTV_Rank, TPTV_initial_rank, TPTV_maxIter, TPTV_lambdas}, ... % TPTV
};

names_params = { ...
    {'epsilon', 'alpha', 'beta', 'boundary', 'blocksize', 'shiftstep', ...
        'maxiter', 'stopcri', 'tmp_save_name'}, ... % S3TTV
    {'epsilon', 'alpha', 'beta', 'boundary', 'maxiter', 'stopcri'}, ... % SSTV
    {'epsilon', 'alpha', 'beta', 'boundary', 'L', 'omega', 'maxiter', 'stopcri'}, ... % HSSTV_L1
    {'epsilon', 'alpha', 'beta', 'boundary', 'L', 'omega', 'maxiter', 'stopcri'}, ... % HSSTV_L12
    {'eps_L2ball', 'alpha', 'beta', 'eps_L10ball', 'stepsize_reduction', ...
        'boundary', 'maxiter', 'stopcri'}, ... % l0l1HTV
    {'epsilon', 'alpha', 'beta', 'boundary', 'blocksize', 'shiftstep', ...
        'maxiter', 'stopcri', 'tmp_save_name'}, ... % STV
    {'epsilon', 'alpha', 'beta', 'boundary', 'blocksize', 'shiftstep', ...
        'maxiter', 'stopcri', 'tmp_save_name'}, ... % SSST
    {'tau', 'lambda', 'rank'}, ... % LRTDTV
    {'Rank', 'initial_rank', 'maxIter', 'lambda'}, ... % TPTV
};


%% Running methods
for idx_method = idc_methods
name_method = names_methods{idx_method};
params_tmp = params_all{idx_method};
name_params = names_params{idx_method};

i_method = i_method + 1;


[params_comb, num_params_comb] = ParamsList2Comb(params_tmp);

for idx_params_comb = 1:num_params_comb

params = struct();
for idx_params = 1:numel(name_params)
    % Assigning parameters to the structure
    params.(name_params{idx_params}) = params_comb{idx_params_comb}{idx_params};
end

switch idx_method
    case 1
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'bl', num2str(params.blocksize(1)), '_st', num2str(params.shiftstep(1)), ...
            '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)); % S3TTV
    case 2
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            '_stop1e-', num2str(stopcri_index)); % SSTV
    case 3
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'o', num2str(params.omega), '_stop1e-', num2str(stopcri_index)); % HSSTV_L1
    case 4
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'o', num2str(params.omega), '_stop1e-', num2str(stopcri_index)); % HSSTV_L12
    case 5
        name_params_savetext = append(...
            'e', num2str(params.eps_L2ball), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'sr', num2str(params.stepsize_reduction), '_th', num2str(l0l1HTV_L10ball_th), '_maxiter', num2str(maxiter)); % l0l1HTV
    case 6
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'bl', num2str(params.blocksize(1)), '_st', num2str(params.shiftstep(1)), ...
            '_stop1e-', num2str(stopcri_index)); % STV
    case 7
        name_params_savetext = append(...
            'e', num2str(params.epsilon), '_a', num2str(params.alpha), '_b', num2str(beta), ...
            'bl', num2str(params.blocksize(1)), '_st', num2str(params.shiftstep(1)), ...
            '_stop1e-', num2str(stopcri_index)); % SSST
    case 8
        name_params_savetext = append('l', num2str(params.lambda), ...
            '_r', num2str(params.rank(1)), '_stop1e-', num2str(stopcri_index)); % LRTDTV
    case 9
        name_params_savetext = append('maxiter', num2str(params.maxIter), ...
            '_l', num2str(params.lambda)); % TPTV
end



fprintf('\n~~~ SETTINGS ~~~\n');
fprintf('Method: %s\n', name_method);
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);
fprintf('Parameter settings: %s\n', name_params_savetext)
fprintf('Methods: (%d/%d), Cases: (%d/%d), Params:(%d/%d)\n', ...
    i_method, numel(idc_methods), idx_exp, total_exp, idx_params_comb, num_params_comb);

[HSI_restored, removed_noise, other_result]...
    = func_methods{idx_method}(params);


% Plotting results
mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
mssim  = calc_MSSIM(HSI_restored, HSI_clean);
ergas  = calc_ERGAS(HSI_restored, 1); % GSD ratio = 1 for recovery problem;
sam    = calc_SAM(HSI_restored, HSI_clean);

fprintf('~~~ RESULTS ~~~\n');
fprintf('MPSNR: %#.4g\n', mpsnr);
fprintf('MSSIM: %#.4g\n', mssim);
fprintf('ERGAS: %#.4g\n', ergas);
fprintf('SAM  : %#.4g\n', sam);

[psnr_per_band, ssim_per_band] = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


% Saving each result
save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

mkdir(save_folder_name);

save(append(save_folder_name, 'image_result.mat'), ...
    'HSI_clean', 'HSI_noisy', 'hsi', 'image', ...
    'HSI_restored', 'removed_noise', ...
    '-v7.3', '-nocompression' ...
);


save(append(save_folder_name, 'other_result.mat'), ...
    'params', 'other_result', ...
    '-v7.3', '-nocompression' ...
);

close all

end

end
end

fprintf('******* finis *******\n');