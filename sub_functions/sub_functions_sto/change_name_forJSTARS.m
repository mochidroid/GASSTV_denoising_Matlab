clear;
close all;

addpath(genpath('sub_functions'))


%% Generating observation
noise_conditions = { ...
    {0.05, 0.05, 0, 0}, ... % g0.05 ps0.05 tint0
    {0.1, 0.05, 0, 0}, ... % g0.1 ps0.05 tint0
    {0.05, 0, 0.05, 0.5}, ... % g0.05 ps0 tint0.5
    {0.1, 0, 0.05, 0.5}, ... % g0.1 ps0 tint0.5
    {0.05, 0.05, 0.05, 0.5}, ... % g0.05 ps0.05 tint0.5
    {0.1, 0.05, 0.05, 0.5}, ... % g0.1 ps0.05 tint0.5
};

num_noise_conditions = size(noise_conditions, 2);

images = {...
    'JasperRidge', ...
    'PaviaU120' ...
};

num_images = numel(images);


for idx_noise_condition = 1:num_noise_conditions
for idx_image = 1:num_images
fprintf('******* initium *******\n');
deg_tmp.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg_tmp.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg_tmp.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg_tmp.stripe_intensity    = noise_conditions{idx_noise_condition}{4};

image = images{idx_image};

[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg_tmp] = Generate_obsv_for_denoising(HSI_clean, deg_tmp);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

%% Setting common parameters
% rho = 0.9;
rho = 0.95;
% rho = 0.95;



% boundary = 'Neumann';
boundary = 'circulant';

blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

stopcri_index = 5;

maxiter = 20000;
% maxiter = 10;
% maxiter = 1;


%% Setting unique parameters
% HSSTV
HSSTV_L1_L = 'L1';
HSSTV_L12_L = 'L12';
HSSTV_omega = 0.05;


% l0l1HTV
l0l1HTV_stepsize_reduction = 0.999;
% l0l1HTV_stepsize_reduction = 0.9999;
% l0l1HTV_stepsize_reduction = 1;



% % LRTDTV
% LRTDTV_tau = 1;
% LRTDTV_lambda = 100;
% LRTDTV_rank = [51,51,10];


% % TPTV
% TPTV_Rank = [7,7,5];
% TPTV_initial_rank = 2;
% TPTV_maxIter = 50;
% TPTV_lambda = 4e-3*sqrt(hsi.n1*hsi.n2);


%% Selecting methods
% name_method = 'SSTV';
% name_method = 'HSSTV_L12';
% name_method = 'l0l1HTV';
% name_method = 'STV';
% name_method = 'SSST';
% name_method = 'S3TTV';
% name_method = 'LRTDTV';
name_method = 'TPTV';


% name_params_savetext = append('r', num2str(rho), '_stop1e-', num2str(stopcri_index));
% name_params_savetext = append('o', num2str(HSSTV_omega), '_r', num2str(rho), ...
%         '_stop1e-', num2str(stopcri_index)); % HSSTV_L1 & L12
% name_params_savetext = append('sr', num2str(l0l1HTV_stepsize_reduction), ...
%         '_r', num2str(rho), '_maxiter', num2str(maxiter)); % l0l1HTV
% name_params_savetext = append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
%         '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)); % STV & SSST & S3TTV
% name_params_save_text = append('stop1e-', num2str(stopcri_index)); % LRTDTV
name_params_save_text = append(''); % TPTV


fprintf('~~~ SETTINGS ~~~\n');
fprintf('Method: %s\n', name_method);
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);
fprintf('Gaussian sigma: %g\n', deg_tmp.gaussian_sigma);
fprintf('Sparse rate: %g\n', deg_tmp.sparse_rate);
fprintf('Stripe rate: %g\n', deg_tmp.stripe_rate);
fprintf('Stripe intensity: %g\n', deg_tmp.stripe_intensity);

%% Loading prev names
deg_prev = deg_tmp;
if deg_prev.sparse_rate == 0
    deg_prev.is_sparse = 0;
else
    deg_prev.is_sparse = 1;
end

if deg_prev.stripe_intensity == 0
    deg_prev.is_stripe = 0;
else
    deg_prev.is_stripe = 1;
end

save_prev_folder_name = Output_save_prev_name(deg_prev, image);

% save_prev_file_name = save_prev_folder_name + "/" + name_method + ...
%     "_cir_PsPPDS_" + image + ...
%     "_r" + num2str(rho) + ...
%     "_stop1e-" + num2str(stopcri_index) + ".mat";

% save_prev_file_name = save_prev_folder_name + "/" + name_method + ...
%     "_cir_PsPPDS_" + image + ...
%     "_sr" + num2str(l0l1HTV_stepsize_reduction) + ...
%     "_r" + num2str(rho) + ...
%     "_maxiter" + num2str(maxiter) + ".mat";

% save_prev_file_name = save_prev_folder_name + "/" + name_method + ...
%     "_cir_PsPPDS_" + image + ...
%     "_bl" + num2str(blocksize(1)) + "_st" + num2str(shiftstep(1)) + ...
%     "_r" + num2str(rho) + ...
%     "_stop1e-" + num2str(stopcri_index) + ".mat"; % STV & SSST & S3TTV

% save_prev_file_name = save_prev_folder_name + "/" + name_method + "_" + ...
%     image + "_stop1e-" + num2str(stopcri_index) + ".mat";

save_prev_file_name = save_prev_folder_name + "/" + name_method + "_" + ...
    image + ".mat";

load(save_prev_file_name)

clear deg
deg = deg_tmp;


%% Modifing prev names
HSI_restored = TPTV_restored_HSI;
% removed_noise = LRTDTV_removed_noise;
% removed_noise.all_noise = S;
% other_result = LRTDTV_per_iter;
% other_result.iteration = LRTDTV_iteration;
% other_result.time = LRTDTV_time;
% other_result.out_value = out_value;

removed_noise.all_noise = HSI_noisy - HSI_restored;



%% Plotting results
mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
mssim  = calc_MSSIM(HSI_restored, HSI_clean);
ergas  = calc_ERGAS(HSI_restored, HSI_clean, 1); % GSD ratio = 1 for recovery problem;
sam    = calc_SAM(HSI_restored, HSI_clean);

[psnr_per_band, ssim_per_band] = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


%% Saving results
save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_tint', num2str(deg.stripe_intensity), '/', ...
    name_method, '/', ...
    name_params_save_text, '/' ...   
);

mkdir(save_folder_name);

save(append(save_folder_name, 'image_result.mat'), ...
    'HSI_clean', 'HSI_noisy', 'hsi', 'deg', 'image', ...
    'HSI_restored', 'removed_noise', ...
    '-v7.3', '-nocompression' ...
);

save(append(save_folder_name, 'metric_vals.mat'), ...
    'mpsnr', 'mssim', 'ergas', 'sam', ...
    'psnr_per_band', 'ssim_per_band', ...
    '-v7.3', '-nocompression' ...
);

% save(append(save_folder_name, 'other_result.mat'), ...
%     'params', 'other_result', ...
%     '-v7.3', '-nocompression' ...
% );

save(append(save_folder_name, 'other_result.mat'), ...
    'params', ...
    '-v7.3', '-nocompression' ...
);

fprintf('******* finis *******\n');


end
end
