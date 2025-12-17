clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Switching operator
is_show_cropped_image = 1;
is_show_HSI = 1;  
% is_output_image = 1;


%% Setting common parameters
images = {...
    "IndianPines", ...
    "Suwannee", ...
};

% Radius parameters; 1st: IndianPines, 2nd: Suwannee
epsilon = {30, 100};
alpha = {200, 800};
beta = {100, 5000};

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;


load("dir_save_comp_folder.mat", "dir_save_comp_folder");



%% Setting output parameters
idx_output = 1;

switch idx_output
    case 1
        idx_image = 1;

        % gain_restoration = 1;
        gain_restoration = 1.5; 
        
        gain_diff = 5;

        save_band = 32; % IndianPines
        % save_band = 40; % IndianPines
        % save_band = 88; % IndianPines

        crop_start_pos = [2, 30];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';
        
        % arrow_head_pos = [20, 8];
        % % arrow_head_pos = [38, 74];
        % % arrow_head_pos = [53, 64];
        % arrow_length = 15;
        % arrow_handle_width = 3;
        % arrow_head_width = 3;
        % arrow_dir_tblr = "l";
        % arrow_methods_idc = [];

    case 2
        idx_image = 2;

        gain_restoration = 1;
        % gain_restoration = 2;
        
        gain_diff = 8;

        save_band = 196; % Suwannee
        
        % crop_start_pos = [27, 80];
        crop_start_pos = [30, 80];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'bl';
        
        arrow_head_pos = [20, 8];
        % arrow_head_pos = [38, 74];
        % arrow_head_pos = [53, 64];
        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_dir_tblr = "l";
        arrow_methods_idc = [];

end

 
%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "output_name", "SSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% l0l1HTV
l0l1HTV.stepsize_reduction = 0.999; % 0.999
l0l1HTV.L10ball_th = 0.03; % 0.02, 0.03, 0.04

methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "output_name", "$\llHTV$", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_sr%.5g_th%.2f_maxiter%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            l0l1HTV.stepsize_reduction, l0l1HTV.L10ball_th, maxiter), ...
    "line_style", "--", ...
    "enable", false ...
);

% HSSTV_L1
HSSTV.omega = 0.05;

methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "output_name", "HSSTV1", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "output_name", "HSSTV2", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% TPTV
TPTV.Rank = [7,7,5];
TPTV.initial_rank = 2;
TPTV.maxIter = 50; % 50, 100
TPTV.lambdas = 1e-3;% 5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2

methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "output_name", "TPTV", ...
    "get_params_savetext", ...
        sprintf("maxiter%d_l%.4g", TPTV.maxIter, TPTV.lambdas), ...
    "line_style", "--", ...
    "enable", true ...
);

% QRNN3D
QRNN3D.ckpt = "complex"; % "complex", "paviaft"

methods_info(end+1) = struct( ...
    "name", "QRNN3D", ...
    "output_name", "QRNN3D", ...
    "get_params_savetext", ...
        sprintf("p_%s_0normalized", QRNN3D.ckpt), ...
    "line_style", "--", ...
    "enable", true ...
);


% GASSTV_CondatVu
GASSTV_CondatVu.sigma_sp = "90"; % {"90", "med"}
GASSTV_CondatVu.sigma_l = "med";
GASSTV_CondatVu.lambda_rho_sp = 0.9; %{0.9, 1, 1.1}
GASSTV_CondatVu.lambda2 = 0.1; % {0.1, 1, 10};
GASSTV_CondatVu.num_segments = 4;
GASSTV_CondatVu.order_filt = 5; % {5, 1};
% GASSTV_CondatVu.order_filt = {5};
methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "output_name", "GASSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_l%.2g_%.2g_sig%s_%s_ns%d_fil%d_stop1e-%d", ...
        epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
        GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, ...
        GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, ...
        stopcri_idx), ...
    "line_style", "-", ...
    "enable", true ...
);

% methods_info(end+1) = struct( ...
%     "name", "GeoSSTV", ...
%     "output_name", "GeoSSTV_v2", ...
%     "get_params_savetext", ...
%         sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
%             110, 850, beta{idx_image}, ...
%             GeoSSTV.lambda, stopcri_idx), ...
%     "line_style", "-", ...
%     "enable", true ...
% );

methods_info = methods_info([methods_info.enable]); % removing true methods
num_methods = numel(methods_info);


%% Loading clean and noisy HS images
image = images{idx_image};

[HSI_noisy, hsi] = Load_real_HSI(image);
HSI_noisy = single(HSI_noisy);


fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);


%% Cropping clean and noisy HS images
image_noisy = HSI_noisy(:,:,save_band)*gain_restoration;
[image_noisy, crop_image_noisy] = Crop_Sep_image(image_noisy, crop_start_pos, crop_size);


% if ~isempty(arrow_methods_idc)
%     image_noisy = Embed_arrow(image_noisy, ...
%                 arrow_head_pos, arrow_length, ...
%                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
% end


%% Loading and cropping result images
for idx_method = 1:num_methods 
    name_method = methods_info(idx_method).name;

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        name_method ...
    );

    params_savetext = methods_info(idx_method).get_params_savetext;

    load(fullfile(dir_result_folder, append(params_savetext, ".mat")), ...
        "HSI_restored");

    methods_info(idx_method).HSI_restored = HSI_restored;


    % Cropping result images
    image_restored = HSI_restored(:,:,save_band)*gain_restoration;
    [image_restored, crop_image_restored] = Crop_Sep_image(image_restored, crop_start_pos, crop_size);

    
    % if find(arrow_methods_idc == idx_method)
    %     image_restored = Embed_arrow(image_restored, ...
    %                 arrow_head_pos, arrow_length, ...
    %                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    % end

    methods_info(idx_method).image_restored         = image_restored;
    methods_info(idx_method).crop_image_restored    = crop_image_restored;

end


%% Showing cropped images
if exist("is_show_cropped_image", "var") && is_show_cropped_image == 1
    cat_restored_image = image_noisy;
    cat_crop_restored_image = crop_image_noisy;

    for idx_method = 1:num_methods 
        cat_restored_image = cat(2, cat_restored_image, methods_info(idx_method).image_restored);
        cat_crop_restored_image = cat(2, cat_crop_restored_image, methods_info(idx_method).crop_image_restored);
    end

    figure;
    imshow(cat_restored_image);

    figure;
    imshow(cat_crop_restored_image);

end

%% Showing HSI
if exist("is_show_HSI", "var") && is_show_HSI == 1
    cat_HSI = HSI_noisy;

    for idx_method = 1:num_methods
        cat_HSI = cat(2, cat_HSI, methods_info(idx_method).HSI_restored);
    end

    cat_diff_GeoSSTV = abs(repmat(methods_info(num_methods).HSI_restored, [1, num_methods+1, 1]) - cat_HSI) * gain_diff;

    implay(cat(1, cat_HSI, cat_diff_GeoSSTV));
end


%% Output restored image
if exist("is_output_image", "var") && is_output_image == 1
    dir_output_root = fullfile(...
        dir_save_comp_folder, ...
        "GeoSSTV_slide_defence", ...
        sprintf("%s_b%d", image, save_band));
    dir_output_result_folder = fullfile(...
        dir_output_root, ...
        "restored_image");
    dir_output_crop_result_folder = fullfile(...
        dir_output_root, ...
        "crop_image");
    mkdir(dir_output_result_folder)
    mkdir(dir_output_crop_result_folder)


    imwrite(image_noisy, ...
        fullfile(dir_output_result_folder, ...
            "image_noisy.png"), ...
        "BitDepth", 8);

    imwrite(crop_image_noisy, ...
        fullfile(dir_output_crop_result_folder, ...
            "crop_image_noisy.png"), ...
        "BitDepth", 8);
    

    for idx_method = 1:num_methods 
        name_method = methods_info(idx_method).name;
        image_restored = methods_info(idx_method).image_restored;
        crop_image_restored = methods_info(idx_method).crop_image_restored;

        imwrite(image_restored, ...
            fullfile(...
                dir_output_result_folder, ...
                sprintf("image_%s.png", name_method)), ...
            "BitDepth", 8);
        imwrite(crop_image_restored, ...
            fullfile(...
                dir_output_crop_result_folder, ...
                sprintf("crop_image_%s.png", name_method)), ...
            "BitDepth", 8);
    end

    fprintf("save dir: %s\n", dir_output_result_folder)
    fprintf("save dir: %s\n", dir_output_crop_result_folder)
end
