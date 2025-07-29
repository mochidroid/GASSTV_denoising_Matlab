clear
close all

%%
% dir_save_ori_folder = fullfile("H:", "マイドライブ", "Matlab_sto", ...
%     "S3TTV_for_JSTARS", "result");
load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder");


% names_methods = { ...
%     "SSTV", ...
%     "HSSTV_L1", ...
%     "HSSTV_L12", ...
%     "l0l1HTV", ...
%     "GASSTV", ...
%     "GASSTV_noisy", ...
%     "GASSTV_Const", ...
% };
% 
% idc_methods = 1:numel(names_methods);
% idc_methods = 1:4;
% idc_methods = 1;


images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

idc_images = 1:numel(images);
% idc_images = 2:3;


noise_conditions = { ...
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0,     0,     0.05,  0.5,  0},     ... % g0 ps0 pt0.05
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = [2:3, 5:8];
% idc_noise_conditions = [7:8];
% idc_noise_conditions = 1;


% 1. フォルダの名前を読み込む
% 2. どこかに格納
% 3. dir_


%%
for idx_image = idc_images
for idx_noise_condition = idc_noise_conditions
clear HSI_clean HSI_noisy ...
    HSI_restored removed_noise ...
    val_mpsnr val_mssim val_sam ...
    vals_psnr_per_band vals_ssim_per_band ...
    params other_result

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

dir_save_result_folder = fullfile(dir_save_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)) ...
);

dir_save_comp_result_folder = fullfile(dir_save_comp_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)) ...
);


fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Image: %s\n", image);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n\n", deg.deadline_rate);


names_methods = dir(dir_save_result_folder);
names_methods = {names_methods.name};
names_methods = names_methods(3:end);


%% 
for idx_method = 1:numel(names_methods)
    name_method = names_methods{idx_method};

    fprintf("%s: ", name_method);

    dir_method_folder = fullfile(dir_save_result_folder, name_method);
    dir_method_comp_folder = fullfile(dir_save_comp_result_folder, name_method);

    mkdir(dir_method_comp_folder)

        
    names_params = dir(fullfile(dir_method_folder, '*.mat'));
    names_params = {names_params.name};

    for idx_params = 1:numel(names_params)
        name_params = names_params(idx_params);

        dir_full_file = fullfile(dir_method_folder, name_params);
        dir_compact_file = fullfile(dir_method_comp_folder, name_params);

        load(dir_full_file, "HSI_restored", "params", "val_mpsnr", "val_mssim", "val_sam");
        
        save(dir_compact_file, "HSI_restored", "params", "val_mpsnr", "val_mssim", "val_sam", ...
            "-v7.3", "-nocompression" ...
            );
    end

    fprintf("Complete\n");

end

end
end


