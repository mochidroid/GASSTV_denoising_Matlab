clear;
close all;

addpath(genpath('sub_functions'))
addpath('func_metrics')


%% Selecting methods
names_methods = { ...
    'SSTV', ...
    'HSSTV_L1', ...
    'HSSTV_L12', ...
    'l0l1HTV', ...
    'STV', ...
    ...
    'SSST', ...
    'LRTDTV', ...
    'TPTV', ...
    'FGSLR', ...
    'LRTDGS', ...
    ...
    'FRCTR_PnP', ...
    'WMLRTR', ...
    'NLSSR', ...
    'FastHyMix', ...
    'S3TTV', ...
};


% idc_methods = 1:numel(names_methods);
idc_methods = [1:9, 14, numel(names_methods)];
% idc_methods = [9, 10, 11];
% idc_methods = 10;


idc_methods_best_param = 9:13;

LineStyles_methods = { ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '--', ...
    '-', ...
};


is_show_HSI = 1;
% is_plot_psnr_and_ssim_per_band = 1;


%% Selecting conditions
noise_conditions = { ...
    % {0.1,   0,      0,      0},     ... % g0.1 ps0 pt0
    % {0,     0,      0.05,   0.5},   ... % g0 ps0 pt0.05
    {0.05,  0.05,   0,      0},     ... % g0.05 ps0.05 pt0
    {0.1,   0.05,   0,      0},     ... % g0.1 ps0.05 pt0
    {0.05,  0,      0.05,   0.5},   ... % g0.05 ps0 pt0.05
    {0.1,   0,      0.05,   0.5},   ... % g0.1 ps0 pt0.05
    {0.05,  0.05,   0.05,   0.5},   ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,   0.05,   0.5},   ... % g0.1 ps0.05 pt0.05
};

idx_noise_condition = 6;
% idc_noise_conditions = 1:size(noise_conditions, 2);

images = {...
    'JasperRidge', ...
    'PaviaU120', ...
    'Beltsville', ...
};

idx_image = 1;
% idc_images = 1:numel(images);
% idc_images = 1:2;

% for idx_image = idc_images
% for idx_noise_condition = idc_noise_conditions
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

%% Setting common parameters
% rho = 0.9;
rho = 0.95;
% rho = 0.95;

% if idx_noise_condition == 1 || idx_noise_condition == 2 
%     rho = 0.98;
% end


% boundary = 'Neumann';
boundary = 'circulant';

blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;

maxiter = 20000;
% maxiter = 10;
% maxiter = 1;


%% Setting unique parameters
% S3TTV
save_S3TTV_Calc_FV_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    'S3TTV/Calc_FV' ...
);

% Check if Calc_FV folder exists
if isfolder(save_S3TTV_Calc_FV_folder_name)
    name_param_save_text_S3TTV = ...
        append('Calc_FV/', ...
            'bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
            '_r', num2str(rho), '_stop1e-', num2str(stopcri_index) ...
        );
else
    name_param_save_text_S3TTV = ...
        append(...
            'bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
            '_r', num2str(rho), '_stop1e-', num2str(stopcri_index) ...
        );
end


% HSSTV
HSSTV_L1_L = 'L1';
HSSTV_L12_L = 'L12';
HSSTV_omega = 0.05;


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


% LRTDTV
LRTDTV_tau = 1;
LRTDTV_lambda = 100;
LRTDTV_rank = [51,51,10];


% TPTV
TPTV_Rank = [7,7,5];
TPTV_initial_rank = 2;
TPTV_maxIter = 50;
TPTV_lambda = 4e-3*sqrt(hsi.n1*hsi.n2);

% GLF
GLF_p_subspace = 8;

% FastHyMix
FastHyMix_k_subspace = 8;


%% Cordinating parameters
names_params_savetext = { ...
    append('r', num2str(rho), '_stop1e-', num2str(stopcri_index)), ... % SSTV
    append('o', num2str(HSSTV_omega), '_r', num2str(rho), ...
        '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L1
    append('o', num2str(HSSTV_omega), '_r', num2str(rho), ...
        '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L12
    append('sr', num2str(l0l1HTV_stepsize_reduction), ...
        '_r', num2str(rho), '_maxiter', num2str(maxiter)) ... % l0l1HTV
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)), ... % STV
    ...
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)), ... % SSST
    append('stop1e-', num2str(stopcri_index)), ... % LRTDTV
    append(''), ... % TPTV
    {}, ... % FGSLR
    {}, ... % LRTDGS
    ...
    {}, ... % FRCTR_PnP
    {}, ... % WMLRTR
    {}, ... % NLSSR
    append('sub', num2str(FastHyMix_k_subspace)), ... % FastHyMix
    name_param_save_text_S3TTV, ... % S3TTV
};

% Required best parameters
for idx_method_best_param = idc_methods_best_param
name_method_best_param = names_methods{idx_method_best_param};

save_best_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method_best_param, '/', ...
    'best_params/' ...   
);

% Get the contents of the directory
contents = dir(save_best_folder_name);
name_param_save_text_best_param = '';

for i = 1:length(contents)
    % Exclude '.' and '..' and check if the item is a folder
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        names_params_savetext{idx_method_best_param} = append('best_params/' , contents(i).name);
        break; % Exit the loop since there is only one folder
    end
end

end

%% Loading restored HS images
cat_HSI = cat(2, HSI_clean, HSI_noisy);
cat_diff = cat(2, zeros(hsi.sizeof), abs(HSI_clean - HSI_noisy));

list_psnr_per_band = [];
list_ssim_per_band = [];

fprintf('~~~ SETTINGS ~~~\n');
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);
fprintf('Gaussian sigma: %g\n', deg.gaussian_sigma);
fprintf('Sparse rate: %g\n', deg.sparse_rate);
fprintf('Stripe rate: %g\n', deg.stripe_rate);
fprintf('Stripe intensity: %g\n', deg.stripe_intensity);


mpsnr_noisy  = calc_MPSNR(HSI_noisy, HSI_clean);
mssim_noisy  = calc_MSSIM(HSI_noisy, HSI_clean);
ergas_noisy  = calc_ERGAS(HSI_noisy, HSI_clean, 1); % GSD ratio = 1 for recovery problem;
sam_noisy    = calc_SAM(HSI_noisy, HSI_clean);

fprintf('~~~ RESULTS ~~~\n');
fprintf('%s  \t MPSNR\t MSSIM\t SAM\n', blanks(max(strlength(names_methods))));

fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append('noisy', blanks(max(strlength(names_methods)) - 5)), ...
    mpsnr_noisy, mssim_noisy, sam_noisy);


for idx_method = idc_methods
name_method = names_methods{idx_method};
name_params_savetext = names_params_savetext{idx_method};

save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_restored' ...
);

load(append(save_folder_name, 'metric_vals.mat'), ...
    'mpsnr', 'mssim', 'sam', 'psnr_per_band', 'ssim_per_band' ...
);


cat_HSI = cat(2, cat_HSI, HSI_restored);
cat_diff = cat(2, cat_diff, abs(HSI_clean - HSI_restored));

list_psnr_per_band = cat(1, list_psnr_per_band, psnr_per_band);
list_ssim_per_band = cat(1, list_ssim_per_band, ssim_per_band);

fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append(name_method, blanks(max(strlength(names_methods)) - strlength(name_method))), ...
    mpsnr, mssim, sam);


end

%% Showing HSI
if exist('is_show_HSI', 'var') && is_show_HSI == 1
cat_diff = cat_diff * diff_magnification;
show_HSI = cat(1, cat_HSI, cat_diff);
implay(show_HSI)
end


%% Plotting psnr and ssim per band
if exist('is_plot_psnr_and_ssim_per_band', 'var') && is_plot_psnr_and_ssim_per_band == 1
label_style = {'FontSize', 15};


figure(1)
hold on
for i = 1:numel(idc_methods)
    name_method = names_methods{idc_methods(i)};
    LineStyle_method = LineStyles_methods{idc_methods(i)};
    plot_style = {'LineWidth', 2.0, 'LineStyle', LineStyle_method, 'DisplayName', name_method};
    plot(list_psnr_per_band(i,:), plot_style{:})
end
hold off
title('psnr per band')
xlabel('band', label_style{:})
ylabel('psnr', label_style{:})

figure(2)
hold on
for i = 1:numel(idc_methods)
    name_method = names_methods{idc_methods(i)};
    LineStyle_method = LineStyles_methods{idc_methods(i)};
    plot_style = {'LineWidth', 2.0, 'LineStyle', LineStyle_method, 'DisplayName', name_method};
    plot(list_ssim_per_band(i,:), plot_style{:})
end
hold off
title('ssim per band')
xlabel('band', label_style{:})
ylabel('psnr', label_style{:})

end

% end
% end