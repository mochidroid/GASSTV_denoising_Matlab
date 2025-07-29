clear
close all

%%
% dir_save_ori_folder = fullfile("H:", "マイドライブ", "Matlab_sto", ...
%     "S3TTV_for_JSTARS", "result");
dir_save_ori_folder = fullfile("H:", "マイドライブ", "MATLAB_Shere", ...
    "RISSTV_for_TGRS_denoising_result");
load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder");


names_methods = { ...
    "SSTV", ...
    "HSSTV_L1", ...
    "HSSTV_L12", ...
    "l0l1HTV", ...
    "GASSTV", ...
    "GASSTV_noisy", ...
    "GASSTV_Const", ...
};

% idc_methods = 1:numel(names_methods);
% idc_methods = 1:4;
idc_methods = 7;


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
% idc_noise_conditions = [2:3, 5:8];
idc_noise_conditions = [7:8];
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

[HSI_clean, hsi] = Load_HSI(image);
noise_seed = "default";
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

dir_save_ori_result = fullfile(dir_save_ori_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)) ...
);

dir_save_result = fullfile(dir_save_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)) ...
);



for idx_method = idc_methods
    name_method = names_methods{idx_method};
    
    dir_method_folder = fullfile(dir_save_result, name_method);
    % dir_method_comp_folder = fullfile(dir_save_comp_folder, name_method);

    mkdir(dir_method_folder);
    % mkdir(dir_method_comp_folder);
    
    dir_method_ori_folder = fullfile(dir_save_ori_result, name_method);
        
    names_params_tmp = dir(dir_method_ori_folder);
    names_params = {names_params_tmp.name};
    
    names_params = names_params(3:end);
    
    idx = strcmp(names_params, 'best_params') + ...
        strcmp(names_params, 'best_params.mat') + ...
        strcmp(names_params, 'best_params_r0.93.mat') + ...
        strcmp(names_params, 'best_params_r0.95.mat') + ...
        strcmp(names_params, 'best_params_r0.98.mat') + ...
        strcmp(names_params, 'image_result.mat') + ...
        strcmp(names_params, 'metric_vals.mat') + ...
        strcmp(names_params, 'other_result.mat');

    names_params = names_params(~idx);

    for idx_params = 1:numel(names_params)
        name_params = names_params(idx_params);

        load(fullfile(dir_method_ori_folder, name_params, "image_result.mat"))
        load(fullfile(dir_method_ori_folder, name_params, "metric_vals.mat"))
        load(fullfile(dir_method_ori_folder, name_params, "other_result.mat"))



        if exist("mpsnr", "var")
            val_mpsnr = mpsnr;
            val_mssim = mssim;
            val_sam   = sam;
            vals_psnr_per_band = psnr_per_band;
            vals_ssim_per_band = ssim_per_band;
        end

        save(fullfile(dir_method_folder, append(name_params, ".mat")), ...
            "HSI_clean", "HSI_noisy", "hsi", "deg", "image", ...
            "HSI_restored", "removed_noise", ...
            "val_mpsnr", "val_mssim", "val_sam", ...
            "vals_psnr_per_band", "vals_ssim_per_band", ...
            "params", "other_result", ...
            "-v7.3", "-nocompression" ...
        );
    end

end

end
end


