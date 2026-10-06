clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Switching operator
is_show_cropped_image = 1;
% is_show_HSI = 1;
% is_plot_psnr_and_ssim_per_band = 1;
is_output_image = 1;


%% Selecting common parameters
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % g0.05 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
};

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

load("dir_save_comp_folder.mat", "dir_save_comp_folder");

color_types = 'hot';
% color_types = 'gray';

cmap = colormap(color_types);
close


%% Setting output parameters
idx_output = 1;

switch idx_output
    case 1
        idx_image = 1;
        idx_noise_condition = 7;

        % gain_restoration = 1;
        gain_restoration = 1.5;
        
        gain_diff = 5;

        % save_band = 6;
        % save_band = 17;

        save_band = 13;
        % save_band = 35; % Jasper Ridge
        
        % crop_start_pos = [60, 79];
        crop_start_pos = [18, 51];

        % crop_start_pos = [45, 58]; % Honmei

        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';
        
        
        % arrow_head_pos = [10, 90];
        % arrow_dir_tblr = "r";

        % arrow_head_pos = [20, 13];
        % arrow_dir_tblr = "l";

        arrow_head_pos = [40, 16];
        arrow_dir_tblr = "b";

        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_methods_idc = [1:7];

    case 2
        idx_image = 1;
        idx_noise_condition = 7;
        % idx_noise_condition = 4;

        gain_restoration = 1.5;
        % gain_restoration = 2;
        
        gain_diff = 5;

        save_band = 20; % Jasper Ridge
        % save_band = 36; % Jasper Ridge
        % save_band = 49; % Jasper Ridge

        % crop_start_pos = [51, 3];
        % crop_start_pos = [45, 3];

        % crop_start_pos = [74, 12];
        crop_start_pos = [78, 16];
        % crop_start_pos = [110, 32];

        crop_size = [20, 20];
        % crop_size = [25, 25];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';
        
        arrow_head_pos = [20, 8];
        % arrow_head_pos = [38, 74];
        % arrow_head_pos = [53, 64];
        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_dir_tblr = "l";
        arrow_methods_idc = [];

    case 3
        idx_image = 1;
        idx_noise_condition = 7;
        % idx_noise_condition = 4;

        gain_restoration = 1.5;
        % gain_restoration = 2;
        
        gain_diff = 5;

        save_band = 20; % Jasper Ridge
        % save_band = 36; % Jasper Ridge
        % save_band = 49; % Jasper Ridge

        % crop_start_pos = [51, 3];
        % crop_start_pos = [45, 3];

        % crop_start_pos = [74, 12];
        crop_start_pos = [78, 16];
        % crop_start_pos = [110, 32];

        crop_size = [20, 20];
        % crop_size = [25, 25];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';
        
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
    "line_style", "--", ...
    "enable", true ...
);

% l0l1HTV
methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "output_name", "$\llHTV$", ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L1
methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "output_name", "HSSTV1", ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "output_name", "HSSTV2", ...
    "line_style", "--", ...
    "enable", true ...
);

% TPTV
methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "output_name", "TPTV", ...
    "line_style", "--", ...
    "enable", true ...
);

% FastHyMix
methods_info(end+1) = struct( ...
    "name", "FastHyMix", ...
    "output_name", "FastHyMix", ...
    "line_style", "--", ...
    "enable", true ...
);

% QRNN3D
methods_info(end+1) = struct( ...
    "name", "QRNN3D", ...
    "output_name", "QRNN3D", ...
    "line_style", "--", ...
    "enable", false ...
);

% GASSTV_CondatVu
methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "output_name", "GASSTV", ...
    "line_style", "-", ...
    "enable", true ...
);


methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Loading clean and noisy HS images
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


% Calculation noisy results
val_mpsnr_noisy  = calc_MPSNR(HSI_noisy, HSI_clean);
val_mssim_noisy  = calc_MSSIM(HSI_noisy, HSI_clean);


fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Case %g (g%g, ps%g, pt%g, pd%g)\n", ...
    idx_noise_condition, deg.gaussian_sigma, deg.sparse_rate, deg.stripe_rate, deg.deadline_rate);


%% Cropping clean and noisy HS images
image_clean = HSI_clean(:,:,save_band)*gain_restoration;
image_noisy = HSI_noisy(:,:,save_band)*gain_restoration;

[image_clean, crop_image_clean] = Crop_Sep_image(image_clean, crop_start_pos, crop_size);
[image_noisy, crop_image_noisy] = Crop_Sep_image(image_noisy, crop_start_pos, crop_size);

diff_image_noisy_gray = abs(image_noisy - image_clean);
diff_image_noisy = Gray2RGB(diff_image_noisy_gray, cmap);
% S
% if ~isempty(arrow_methods_idc)
%     image_clean = Embed_arrow(image_clean, ...
%                 arrow_head_pos, arrow_length, ...
%                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
%     image_noisy = Embed_arrow(image_noisy, ...
%                 arrow_head_pos, arrow_length, ...
%                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
%     diff_image_noisy = Embed_arrow(diff_image_noisy, ...
%                 arrow_head_pos, arrow_length, ...
%                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
% end


%% Loading and cropping result images
for idx_method = 1:num_methods 
    output_name_method = methods_info(idx_method).name;

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
        output_name_method ...
    );

    load(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
    methods_info(idx_method).get_params_savetext = best_params_savetext;

    load(fullfile(dir_result_folder, append(best_params_savetext, ".mat")), ...
        "HSI_restored", "val_mpsnr", "val_mssim", "val_sam");

    methods_info(idx_method).HSI_restored = HSI_restored;
    methods_info(idx_method).val_mpsnr = val_mpsnr;
    methods_info(idx_method).val_mssim = val_mssim;


    % Cropping result images
    image_restored = HSI_restored(:,:,save_band)*gain_restoration;
    
    diff_image_restored_gray = abs(image_restored - image_clean)*gain_diff;
    diff_image_restored = Gray2RGB(diff_image_restored_gray, cmap);

    [image_restored, crop_image_restored] = Crop_Sep_image(image_restored, crop_start_pos, crop_size);
    
    
    % if find(arrow_methods_idc == idx_method)
    %     image_restored = Embed_arrow(image_restored, ...
    %                 arrow_head_pos, arrow_length, ...
    %                 arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    %     diff_image_restored = Embed_arrow(diff_image_restored, ...
    %             arrow_head_pos, arrow_length, ...
    %             arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    % end


    methods_info(idx_method).image_restored         = image_restored;
    methods_info(idx_method).crop_image_restored    = crop_image_restored;
    methods_info(idx_method).diff_image_restored    = diff_image_restored;
end


%% Plotting metrics
% Initialization
name_length_max = max(strlength([methods_info.name]));
name_params_length_max = max(strlength([methods_info.get_params_savetext]));

fprintf("~~~ RESULTS ~~~\n");
fprintf("%s  \t MPSNR\t MSSIM\n", blanks(name_length_max + name_params_length_max + 2));

% Plotting noisy results
fprintf("%s: \t %#.4g\t %#.4g\n", ...
    append("noisy", blanks(name_length_max + name_params_length_max - 3)), ...
    val_mpsnr_noisy, val_mssim_noisy);

% Plotting restored results
for idx_method = 1:num_methods
    output_name_method         = methods_info(idx_method).name;
    params_save_text    = methods_info(idx_method).get_params_savetext;
    val_mpsnr           = methods_info(idx_method).val_mpsnr;
    val_mssim           = methods_info(idx_method).val_mssim;

    fprintf("%s(%s): \t %#.4g\t %#.4g\n", ...
        append(output_name_method, blanks(name_length_max - strlength(output_name_method))), ...
        append(params_save_text, blanks(name_params_length_max - strlength(params_save_text))), ...
        val_mpsnr, val_mssim);
end

fprintf("\n\n")


%% Showing cropped images
if exist("is_show_cropped_image", "var") && is_show_cropped_image == 1
    cat_restored_image = cat(2, image_clean, image_noisy);
    cat_crop_restored_image = cat(2, crop_image_clean, crop_image_noisy);
    cat_diff_image = cat(2, zeros([hsi.n1, hsi.n2, 3]), diff_image_noisy);

    for idx_method = 1:num_methods 
        cat_restored_image = cat(2, cat_restored_image, methods_info(idx_method).image_restored);
        cat_crop_restored_image = cat(2, cat_crop_restored_image, methods_info(idx_method).crop_image_restored);
        cat_diff_image = cat(2, cat_diff_image, methods_info(idx_method).diff_image_restored);
    end

    figure;
    imshow(cat_restored_image);

    figure;
    imshow(cat_crop_restored_image);
    
    figure;
    imshow(cat_diff_image)
end


%% Showing HSI
if exist("is_show_HSI", "var") && is_show_HSI == 1
    cat_HSI = cat(2, HSI_clean, HSI_noisy);

    for idx_method = 1:num_methods
        cat_HSI = cat(2, cat_HSI, methods_info(idx_method).HSI_restored);
    end

    cat_diff = abs(repmat(HSI_clean, [1, num_methods+2, 1]) - cat_HSI) * gain_diff;
    diff_GeoSSTV = abs(HSI_clean - methods_info(num_methods).HSI_restored) * gain_diff;
    cat_diff_GeoSSTV = repmat(diff_GeoSSTV, [1, num_methods+2, 1]) - cat_diff + 0.5;

    implay(cat(1, cat_HSI, cat_diff, cat_diff_GeoSSTV));
end


%% Plotting psnr and ssim per band
if exist("is_plot_psnr_and_ssim_per_band", "var") && is_plot_psnr_and_ssim_per_band == 1
    % Calculating psnr and ssim per band
    vals_psnr_per_band = zeros(num_methods, hsi.n3);
    vals_ssim_per_band = zeros(num_methods, hsi.n3);
    
    for idx_method = 1:num_methods
        [psnr_per_band, ssim_per_band] = ...
            calc_PSNR_SSIM_per_band(methods_info(idx_method).HSI_restored, HSI_clean);
        vals_psnr_per_band(idx_method, :) = psnr_per_band;
        vals_ssim_per_band(idx_method, :) = ssim_per_band;
    end
    
    % Plotting results
    label_style = {"FontSize", 15};
    
    figure()
    hold on
    for idx_method = 1:num_methods
        output_name_method = methods_info(idx_method).name;
        lineStyle_method = methods_info(idx_method).line_style;
        plot_style = {"LineWidth", 2.0, "LineStyle", lineStyle_method, "DisplayName", output_name_method};
        plot(vals_psnr_per_band(idx_method,:), plot_style{:})
    end
    hold off
    title("psnr per band")
    xlabel("band", label_style{:})
    ylabel("psnr", label_style{:})
    
    figure()
    hold on
    for idx_method = 1:num_methods
        output_name_method = methods_info(idx_method).name;
        lineStyle_method = methods_info(idx_method).line_style;
        plot_style = {"LineWidth", 2.0, "LineStyle", lineStyle_method, "DisplayName", output_name_method};
        plot(vals_ssim_per_band(idx_method,:), plot_style{:})
    end
    hold off
    title("ssim per band")
    xlabel("band", label_style{:})
    ylabel("ssim", label_style{:})

end

%% Output restored image
if exist("is_output_image", "var") && is_output_image == 1
    dir_output_root = fullfile(...
        dir_save_comp_folder, ...
        "GASSTV_slide_defence", ...
        sprintf("%s_case%d_b%d", image, idx_noise_condition, save_band));
    dir_output_result_folder = fullfile(...
        dir_output_root, ...
        "restored_image");
    dir_output_crop_result_folder = fullfile(...
        dir_output_root, ...
        "crop_image");
    % dir_output_diff_folder = fullfile(...
    %     dir_output_root, ...
    %     sprintf("diff_image_%s", color_types));
    mkdir(dir_output_result_folder)
    mkdir(dir_output_crop_result_folder)
    % mkdir(dir_output_diff_folder)

    % Saving clean and noisy images
    imwrite(image_clean, ...
        fullfile(dir_output_result_folder, ...
            "image_clean.png"), ...
        "BitDepth", 8);

    imwrite(image_noisy, ...
        fullfile(dir_output_result_folder, ...
            "image_noisy.png"), ...
        "BitDepth", 8);

    imwrite(crop_image_clean, ...
        fullfile(dir_output_crop_result_folder, ...
            "crop_image_clean.png"), ...
        "BitDepth", 8);

    imwrite(crop_image_noisy, ...
        fullfile(dir_output_crop_result_folder, ...
            "crop_image_noisy.png"), ...
        "BitDepth", 8);
    

    for idx_method = 1:num_methods 
        output_name_method = methods_info(idx_method).output_name;
        image_restored = methods_info(idx_method).image_restored;
        crop_image_restored = methods_info(idx_method).crop_image_restored;
        % diff_image_restored = methods_info(idx_method).diff_image_restored;

        imwrite(image_restored, ...
            fullfile(...
                dir_output_result_folder, ...
                sprintf("image_%s.png", output_name_method)), ...
            "BitDepth", 8);
        imwrite(crop_image_restored, ...
            fullfile(...
                dir_output_crop_result_folder, ...
                sprintf("crop_image_%s.png", output_name_method)), ...
            "BitDepth", 8);
        % imwrite(diff_image_restored, ...
        %     fullfile(...
        %         dir_output_diff_folder, ...
        %         sprintf("image_%s.png", name_method)), ...
        %     "BitDepth", 8);
    end

    fprintf("save dir: %s\n", dir_output_result_folder)
    fprintf("save dir: %s\n", dir_output_crop_result_folder)
    % fprintf("save dir: %s\n", dir_output_diff_folder)
end