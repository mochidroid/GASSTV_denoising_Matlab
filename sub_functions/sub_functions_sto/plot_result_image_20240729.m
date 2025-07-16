clear;
close all;

addpath(genpath('sub_functions'))



%% Selecting methods
names_methods = { ...
    'SSTV', ...
    'HSSTV_L1', ...
    'HSSTV_L12', ...
    'l0l1HTV', ...
    'STV', ...
    'SSST', ...
    'LRTDTV', ...
    'TPTV', ...
    'S3TTV', ...
};

% idc_methods = 1:num_methods;
idc_methods = [1:4,9];

num_methods = numel(idc_methods);

%% Selecting noise condition
deg.gaussian_sigma      = 0.1; % standard derivation of Gaussian noise
deg.sparse_rate         = 0.05;
deg.stripe_rate         = 0.05;
deg.stripe_intensity    = 0.5;

% image = 'JasperRidge';
% image = 'PaviaU120';
% image = 'WashingtonDC';
image = 'Salinas';


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


diff_magnification = 7;

%% Setting common parameters
% rho = 0.9;
rho = 0.95;
% rho = 0.95;

epsilon = rho * deg.gaussian_sigma * sqrt(hsi.N);


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
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)), ... % SSST
    append('stop1e-', num2str(stopcri_index)), ... % LRTDTV
    append(''), ... % TPTV
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)) ... % S3TTV
};

%% Loading restored HS images
cat_HSI = cat(2, HSI_clean, HSI_noisy);
cat_diff = cat(2, zeros(hsi.sizeof), abs(HSI_clean - HSI_noisy));

fprintf('~~~ RESULTS ~~~\n');
fprintf('%s  \t MPSNR\t MSSIM\t SAM\n', blanks(max(strlength(names_methods))));

for idx_method = idc_methods
name_method = names_methods{idx_method};
name_params_savetext = names_params_savetext{idx_method};

save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_tint', num2str(deg.stripe_intensity), '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_restored' ...
);

load(append(save_folder_name, 'metric_vals.mat'), ...
    'mpsnr', 'mssim', 'sam' ...
);


cat_HSI = cat(2, cat_HSI, HSI_restored);
cat_diff = cat(2, cat_diff, abs(HSI_clean - HSI_restored));

fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append(name_method, blanks(max(strlength(names_methods)) - strlength(name_method))), ...
    mpsnr, mssim, sam);


end

cat_diff = cat_diff * diff_magnification;
show_HSI = cat(1, cat_HSI, cat_diff);
implay(show_HSI)
