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


is_method = [1, 2, 6, 9];
is_method = [9];
% is_method = 1:numel(names_methods);
num_methods = numel(is_method);

%% Selecting noise condition
deg.gaussian_sigma      = 0.1; % standard derivation of Gaussian noise
deg.sparse_rate         = 0.05;
deg.stripe_rate         = 0.05;
deg.stripe_intensity    = 0.5;

image = 'JasperRidge';
% image = 'PaviaU120';
% image = 'WashingtonDC';


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


fig_font_size = 10;

%% Setting difference operator
D       = @(z) cat(4, z([2:end, 1],:,:) - z, z(:,[2:end, 1],:) - z);
Dt      = @(z) z([end,1:end-1],:,:,1) - z(:,:,:,1) + z(:,[end,1:end-1],:,2) - z(:,:,:,2);
% Dv      = @(z) z([2:end, 1],:,:) - z;
% Dvt     = @(z) z([end,1:end-1],:,:) - z(:,:,:);
Ds      = @(z) z(:, :, [2:end, 1], :) - z;
Dst     = @(z) z(:,:,[end,1:end-1],:) - z(:,:,:,:);


%% Setting common parameters
rho = 0.95;

epsilon = rho * deg.gaussian_sigma * sqrt(hsi.N * (1 - deg.sparse_rate));
alpha = rho * (0.5 * hsi.N * deg.sparse_rate);
beta = rho * hsi.N * (1 - deg.sparse_rate) ...
    * deg.stripe_rate * deg.stripe_intensity / 2;


% boundary = 'Neumann';
boundary = 'circulant';

blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;

maxiter = 20000;


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
DDsU_true   = D(Ds(HSI_clean));
DDsV        = D(Ds(HSI_noisy));

DDsU_true_vec   = abs(DDsU_true(:));
DDsV_vec        = abs(DDsV(:));


figure(1)
    subplot(1, num_methods+2, 1)
    plot(DDsU_true_vec)
    title('Ground-truth')
    ylim([0,1])
    xlabel('Element', 'FontSize', fig_font_size)
    ylabel('Absolute value', 'FontSize', fig_font_size)
    
    % figure(2)
    subplot(1, num_methods+2, 2)
    plot(DDsV_vec)
    title('Observation')
    xlabel('Element', 'FontSize', fig_font_size)
    ylabel('Absolute value', 'FontSize', fig_font_size)
    ylim([0,1])


% for idx_method = 1:num_methods
i = 0;
for idx_method = is_method
name_method = names_methods{idx_method};
name_params_savetext = names_params_savetext{idx_method};

i = i + 1;

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

DDsU        = D(Ds(HSI_restored));
DDsU_vec    = abs(DDsU(:));


subplot(1, num_methods+2, i+2)
    plot(DDsU_vec)
    title(name_method)
    xlabel('Element', 'FontSize', fig_font_size)
    ylabel('Absolute value', 'FontSize', fig_font_size)
    ylim([0,1])


end