clear
close all
addpath(genpath('sub_functions'))

rng('default')

%% Setting conditions
deg.gaussian_sigma      = 0.1; % standard derivation of Gaussian noise
deg.sparse_rate         = 0.05;

deg.stripe_rate         = 0.07;
deg.stripe_intensity    = 0.5;

image = 'JasperRidge';
% image = 'PaviaU120';

% multiple = 1.2;
addition = 0.1;

start_band = 10;


[HSI_clean, hsi] = Load_HSI(image);
% HSI_clean = HSI_clean * multiple;
HSI_clean = HSI_clean + addition;
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


gaussian_noise = deg.gaussian_noise + 0.5;
sparse_noise = deg.sparse_noise + 0.5;
stripe_noise = deg.stripe_noise + 0.5;


%% Loading deep restored HS image
name_method = 'FastHyMix';
FastHyMix_k_subspace = 8;
name_params_savetext = append('sub', num2str(FastHyMix_k_subspace));

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

%% Generating cube
rng('default')
% select_band = 48;

cube_clean      = [];
cube_noisy      = [];
cube_gaussian   = [];
cube_sparse     = [];
cube_stripe     = [];
cube_deep   = [];

for i = start_band:hsi.n3
    cube_clean = cat(3, cube_clean, HSI_clean(:,:,i));
    cube_clean = cat(3, cube_clean, HSI_clean(:,:,i));

    cube_noisy = cat(3, cube_noisy, HSI_noisy(:,:,i));
    cube_noisy = cat(3, cube_noisy, HSI_noisy(:,:,i));

    cube_gaussian = cat(3, cube_gaussian, gaussian_noise(:,:,i));
    cube_gaussian = cat(3, cube_gaussian, gaussian_noise(:,:,i));

    cube_sparse = cat(3, cube_sparse, sparse_noise(:,:,i));
    cube_sparse = cat(3, cube_sparse, sparse_noise(:,:,i));

    cube_stripe = cat(3, cube_stripe, stripe_noise(:,:,i));
    cube_stripe = cat(3, cube_stripe, stripe_noise(:,:,i));

    cube_deep = cat(3, cube_deep, HSI_restored(:,:,i));
    cube_deep = cat(3, cube_deep, HSI_restored(:,:,i));
end

cube_n1 = size(cube_clean, 1);
cube_n2 = size(cube_clean, 2);
cube_n3 = size(cube_clean, 3);

%% Generating flat image
cube_front_clean = cube_clean(:,:,1);
cube_front_noisy = cube_noisy(:,:,1);
cube_front_gaussian = cube_gaussian(:,:,1);
cube_front_sparse = cube_sparse(:,:,1);
cube_front_stripe = cube_stripe(:,:,1);
cube_front_deep = cube_deep(:,:,1);


cube_top_clean = rot90(reshape(cube_clean(1,:,:), [cube_n2, cube_n3]));
cube_top_noisy = rot90(reshape(cube_noisy(1,:,:), [cube_n2, cube_n3]));
cube_top_gaussian = rot90(reshape(cube_gaussian(1,:,:), [cube_n2, cube_n3]));
cube_top_sparse = rot90(reshape(cube_sparse(1,:,:), [cube_n2, cube_n3]));
cube_top_stripe = rot90(reshape(cube_stripe(1,:,:), [cube_n2, cube_n3]));
cube_top_deep = rot90(reshape(cube_deep(1,:,:), [cube_n2, cube_n3]));


cube_side_clean = reshape(cube_clean(:,end,:), [cube_n1, cube_n3]);
cube_side_noisy = reshape(cube_noisy(:,end,:), [cube_n1, cube_n3]);
cube_side_gaussian = reshape(cube_gaussian(:,end,:), [cube_n1, cube_n3]);
cube_side_sparse = reshape(cube_sparse(:,end,:), [cube_n1, cube_n3]);
cube_side_stripe = reshape(cube_stripe(:,end,:), [cube_n1, cube_n3]);
cube_side_deep = reshape(cube_deep(:,end,:), [cube_n1, cube_n3]);

%% Showing flat image
is_show = 1;
if exist('is_show', 'var') && is_show
figure;
imshow(cat(2, ...
    cat(1, cube_top_clean, cube_front_clean), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_clean), ...
        ones(cube_n3+cube_n1, 1), ...
    cat(1, cube_top_noisy, cube_front_noisy), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_noisy), ...
        ones(cube_n3+cube_n1, 1), ...
    cat(1, cube_top_gaussian, cube_front_gaussian), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_gaussian), ...
        ones(cube_n3+cube_n1, 1), ...
    cat(1, cube_top_sparse, cube_front_sparse), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_sparse), ...
        ones(cube_n3+cube_n1, 1), ...
    cat(1, cube_top_stripe, cube_front_stripe), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_stripe), ...
    cat(1, cube_top_deep, cube_front_deep), ...
        cat(1, zeros(cube_n3, cube_n3), cube_side_deep) ))
    
end


%% Saving flat image
is_save = 1;
if exist('is_save', 'var') && is_save
    save_folder_name = append(...
        './result_image_for_JSTARS/HSI_obsv_', image);
    % save_folder_name = './result_slide_image/HSI_cube_ver2';
    mkdir(save_folder_name);
    
    imwrite(cube_front_clean, save_folder_name + "/clean_cube_front.png");
    imwrite(cube_top_clean, save_folder_name + "/clean_cube_top.png");
    imwrite(cube_side_clean, save_folder_name + "/clean_cube_side.png");
    
    imwrite(cube_front_noisy, save_folder_name + "/noisy_cube_front.png");
    imwrite(cube_top_noisy, save_folder_name + "/noisy_cube_top.png");
    imwrite(cube_side_noisy, save_folder_name + "/noisy_cube_side.png");
    
    imwrite(cube_front_gaussian, save_folder_name + "/gaussian_cube_front.png");
    imwrite(cube_top_gaussian, save_folder_name + "/gaussian_cube_top.png");
    imwrite(cube_side_gaussian, save_folder_name + "/gaussian_cube_side.png");
    
    imwrite(cube_front_sparse, save_folder_name + "/sparse_cube_front.png");
    imwrite(cube_top_sparse, save_folder_name + "/sparse_cube_top.png");
    imwrite(cube_side_sparse, save_folder_name + "/sparse_cube_side.png");

    imwrite(cube_front_stripe, save_folder_name + "/stripe_cube_front.png");
    imwrite(cube_top_stripe, save_folder_name + "/stripe_cube_top.png");
    imwrite(cube_side_stripe, save_folder_name + "/stripe_cube_side.png");

    imwrite(cube_front_deep, save_folder_name + "/deep_cube_front.png");
    imwrite(cube_top_deep, save_folder_name + "/deep_cube_top.png");
    imwrite(cube_side_deep, save_folder_name + "/deep_cube_side.png");
end

