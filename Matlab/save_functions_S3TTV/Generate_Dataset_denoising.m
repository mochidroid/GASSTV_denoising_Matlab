clear
close all;

% addpath('function_S3TTV')
addpath(genpath('sub_functions'))

fprintf('******* initium *******\n');
rng('default')

%% Selecting conditions
noise_conditions = { ...
    {0.1, 0, 0, 0}, ... % g0.1 ps0 pt0
    {0, 0, 0.05, 0.5}, ... % g0 ps0 pt0.05
    {0.05, 0.05, 0, 0}, ... % g0.05 ps0.05 pt0
    {0.1, 0.05, 0, 0}, ... % g0.1 ps0.05 pt0
    {0.05, 0, 0.05, 0.5}, ... % g0.05 ps0 pt0.05
    {0.1, 0, 0.05, 0.5}, ... % g0.1 ps0 pt0.05
    {0.05, 0.05, 0.05, 0.5}, ... % g0.05 ps0.05 pt0.05
    {0.1, 0.05, 0.05, 0.5}, ... % g0.1 ps0.05 pt0.05
};

idc_noise_conditions = 1:size(noise_conditions, 2);


images = {...
    'JasperRidge', ...
    'PaviaU120', ...
};

idc_images = 1:2;
% idc_images = 1;


for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Generating observation
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
image = images{idx_image};

[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

% HSI_clean = single(HSI_clean);
% HSI_noisy = single(HSI_noisy);


%% Saving each image



% Saving each result
save_folder_name = append(...
    './Deep/Dataset/' , ...
    image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/' ... 
);

mkdir(save_folder_name);

save(append(save_folder_name, 'HSI_clean.mat'), ...
    'HSI_clean' ...
);

save(append(save_folder_name, 'HSI_noisy.mat'), ...
    'HSI_noisy' ...
);

close all

end
end

fprintf('******* finis *******\n');