clear
close all
addpath(genpath('sub_functions'))

% Case6: Jasper Ridge g0.1 ps0.05 tint05

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
    'FGSLR', ...
    'TPTV', ...
    'FastHyMix', ...
    ...
    'S3TTV', ...
};


idc_methods = 1:numel(names_methods);
% idc_methods = [9, 10, 11];
% idc_methods = 10;


idc_methods_best_param = 8:10;


% is_plot = 1;
is_save = 1;


%% Selecting conditions
noise_condition = {0.1, 0.05, 0.05, 0.5}; % g0.05 ps0.05 tint0

image = 'JasperRidge';
% image = 'PaviaU120';
% image = 'Beltsville';


%% Selecting noise condition
deg.gaussian_sigma      = noise_condition{1};
deg.sparse_rate         = noise_condition{2};
deg.stripe_rate         = noise_condition{3};
deg.stripe_intensity    = noise_condition{4};


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

fprintf('~~~ SETTINGS ~~~\n');
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);
fprintf('Gaussian sigma: %g\n', deg.gaussian_sigma);
fprintf('Sparse rate: %g\n', deg.sparse_rate);
fprintf('Stripe rate: %g\n', deg.stripe_rate);
fprintf('Stripe intensity: %g\n', deg.stripe_intensity);


%% Setting save parameters
restoration_magnification = 1;
% restoration_multiple = 1.5;

diff_magnification = 8;

types = 'hot';
% types = 'gray';
% types = 'turbo';
% types = 'jet';
% types = 'parula';

cmap = colormap(types);

% save_band = 59; %Jasper Ridge1
save_band = 76; % Jasper Ridge
% save_band = 129; %Jasper Ridge2
% save_band = 58; % PaviaU g0.05 ps0.05 tint0
% save_band = 131;

% crop_start_pos = [60, 79];
% crop_start_pos = [38, 52];
crop_start_pos = [70, 57];
% crop_start_pos = [72, 55];
crop_size = [20, 20];
crop_expansion_rate = 2;
crop_embed_tblr = 'bl';
% crop_embed_tblr = 'br';

arrow_head_pos = [20, 8];
% arrow_head_pos = [38, 74];
% arrow_head_pos = [53, 64];
arrow_length = 15;
arrow_handle_width = 3;
arrow_head_width = 3;
arrow_dir_tblr = "l";

% arrow_methods_idc = [7,8,9];
% arrow_methods_idc = [7,8];
arrow_methods_idc = [7, 9, 10];


%% Setting common parameters
% rho = 0.9;
rho = 0.95;
% rho = 0.95;


blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

stopcri_index = 5;

maxiter = 20000;


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
    {}, ... % FGSLR
    {}, ... % TPTV
    {}, ... % FastHyMix
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


%% Ploting cropped images
if exist('is_plot', 'var') && is_plot
%% Cropping clean and noisy images
image_clean = HSI_clean(:,:,save_band)*restoration_magnification;
image_clean = Crop_Embed_image(image_clean, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

image_noisy = HSI_noisy(:,:,save_band)*restoration_magnification;
image_noisy = Crop_Embed_image(image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

diff_image_noisy_gray = abs(image_noisy - image_clean);
diff_image_noisy = Gray2RGB(diff_image_noisy_gray, cmap);
diff_image_noisy = Crop_Embed_image(diff_image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

if ~isempty(arrow_methods_idc)
    image_noisy = Embed_arrow(image_noisy, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    diff_image_noisy = Embed_arrow(diff_image_noisy, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
end

cat_restored_image = cat(2, image_clean, image_noisy);
cat_diff_image = cat(2, zeros([hsi.n1, hsi.n2, 3]), diff_image_noisy);


fprintf('~~~ RESULTS ~~~\n');
fprintf('%s  \t MPSNR\t MSSIM\t SAM\n', blanks(max(strlength(names_methods))));

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


image_restored = HSI_restored(:,:,save_band)*restoration_magnification;
image_restored = Crop_Embed_image(image_restored, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

diff_image_restored_gray = abs(image_restored - image_clean)*diff_magnification;
diff_image_restored = Gray2RGB(diff_image_restored_gray, cmap);
diff_image_restored = Crop_Embed_image(diff_image_restored, ...
                crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

if find(arrow_methods_idc == idx_method)
    image_restored = Embed_arrow(image_restored, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    diff_image_restored = Embed_arrow(diff_image_restored, ...
            arrow_head_pos, arrow_length, ...
            arrow_handle_width, arrow_head_width, arrow_dir_tblr);
end

cat_restored_image = cat(2, cat_restored_image, image_restored);
cat_diff_image = cat(2, cat_diff_image, diff_image_restored);


fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append(name_method, blanks(max(strlength(names_methods)) - strlength(name_method))), ...
    mpsnr, mssim, sam);
end

%% Showing images
figure;
imshow(cat_restored_image);

figure;
imshow(cat_diff_image)


end


%% Saving cropped images
if exist('is_save', 'var') && is_save
save_folder_name_cropped_image = append(...
    './result_image_for_JSTARS/', ...
    'Case8_', image ...
);

save_folder_name_restored_image = append(...
    save_folder_name_cropped_image, '/restored_image/', ...
     'b', num2str(save_band), '_m', num2str(restoration_magnification) ...
);

save_folder_name_diff_image = append(...
    save_folder_name_cropped_image, '/diff_image/', ...
     'b', num2str(save_band), '_m', num2str(restoration_magnification), ...
     '_', types ...
);

mkdir(save_folder_name_restored_image);
mkdir(save_folder_name_diff_image)


%% Saving clean and noisy images
image_clean = HSI_clean(:,:,save_band)*restoration_magnification;
image_clean = Crop_Embed_image(image_clean, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

image_noisy = HSI_noisy(:,:,save_band)*restoration_magnification;
image_noisy = Crop_Embed_image(image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);
image_noisy = Embed_arrow(image_noisy, ...
            arrow_head_pos, arrow_length, ...
            arrow_handle_width, arrow_head_width, arrow_dir_tblr);

diff_image_noisy_gray = abs(image_noisy - image_clean);
diff_image_noisy = Gray2RGB(diff_image_noisy_gray, cmap);
diff_image_noisy = Crop_Embed_image(diff_image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);
diff_image_noisy = Embed_arrow(diff_image_noisy, ...
            arrow_head_pos, arrow_length, ...
            arrow_handle_width, arrow_head_width, arrow_dir_tblr);


imwrite(image_clean, save_folder_name_restored_image + "/image_clean.png");
imwrite(image_noisy, save_folder_name_restored_image + "/image_noisy.png");

imwrite(diff_image_noisy, save_folder_name_diff_image + "/diff_image_noisy.png");


%% Saving each images
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


image_restored = HSI_restored(:,:,save_band)*restoration_magnification;
image_restored = Crop_Embed_image(image_restored, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

diff_image_restored_gray = abs(image_restored - image_clean)*diff_magnification;
diff_image_restored = Gray2RGB(diff_image_restored_gray, cmap);
diff_image_restored = Crop_Embed_image(diff_image_restored, ...
                crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

if find(arrow_methods_idc == idx_method)
    image_restored = Embed_arrow(image_restored, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    diff_image_restored = Embed_arrow(diff_image_restored, ...
            arrow_head_pos, arrow_length, ...
            arrow_handle_width, arrow_head_width, arrow_dir_tblr);
end

imwrite(image_restored, ...
    append(save_folder_name_restored_image, ...
    '/image_', name_method, '.png'...
    ));
imwrite(diff_image_restored, ...
    append(save_folder_name_diff_image, ...
    '/diff_image_', name_method, '.png'...
    ));

end
end