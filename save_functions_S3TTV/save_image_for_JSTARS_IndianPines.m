% close all
addpath(genpath('sub_functions'))

% IndianPines

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


is_plot = 1;
% is_save = 1;


%% Selecting noise condition
image = 'IndianPines';

[HSI_noisy, hsi] = Load_real_HSI(image);
HSI_noisy = single(HSI_noisy);

fprintf('~~~ SETTINGS ~~~\n');
fprintf('Image: %s Size: (%d, %d, %d)\n', image, hsi.n1, hsi.n2, hsi.n3);


%% Setting parameters
% restoration_Boost = 1;
restoration_Boost = 1.5;

save_band = 88;

crop_start_pos = [62, 42];

crop_size = [20, 20];
crop_expansion_rate = 2;
% crop_embed_tblr = 'bl';
crop_embed_tblr = 'br';

arrow_head_pos = {[72, 55]};
% arrow_head_pos = [53, 64];
arrow_length = 15;
arrow_handle_width = 3;
arrow_head_width = 3;
arrow_dir_tblr = {"b"};

arrow_methods_idc = {[1, 3, 4, 6, 7, 8, 9, 10, 11]};


%% Setting common parameters
epsilon = 30;
alpha = 200;
beta = 100;

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


%% Cordinating parameters
names_params_savetext = { ...
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_stop1e-', num2str(stopcri_index)), ... % SSTV
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_o', num2str(HSSTV_omega), '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L1
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_o', num2str(HSSTV_omega), '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L12
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_sr', num2str(l0l1HTV_stepsize_reduction), '_th', num2str(l0l1HTV_L10ball_th), ...
        '_maxiter', num2str(maxiter)), ... % l0l1HTV
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_stop1e-', num2str(stopcri_index)), ... % STV
    ...
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_stop1e-', num2str(stopcri_index)), ... % SSST
    append('l', num2str(LRTDTV_lambda), '_r', num2str(LRTDTV_rank(1)), ...
        '_stop1e-', num2str(stopcri_index)), ... % LRTDTV
    {}, ... % FGSLR
    {}, ... % TPTV
    {}, ... % FastHyMix
    ...
    append('e', num2str(epsilon), '_a', num2str(alpha), '_b', num2str(beta), ...
        '_bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_stop1e-', num2str(stopcri_index)), ... % S3TTV
};


%% Cropping clean and noisy images
image_noisy = HSI_noisy(:,:,save_band)*restoration_Boost;

if ~isempty(arrow_methods_idc{1})
    for i = 1:numel(arrow_methods_idc)
        image_noisy = Embed_arrow(image_noisy, ...
                    arrow_head_pos{i}, arrow_length, ...
                    arrow_handle_width, arrow_head_width, arrow_dir_tblr{i});
    end
end

image_noisy = Crop_Embed_image(image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);


%% Cropping restored images
images_restored = zeros([hsi.n1, hsi.n2, numel(names_methods)]);


for idx_method = idc_methods
name_method = names_methods{idx_method};
name_params_savetext = names_params_savetext{idx_method};


% Required best parameters
if find(idc_methods_best_param == idx_method)
    name_method_best_param = names_methods{idx_method};
    
    save_best_folder_name = append(...
        './result/' , ...
        'denoising_', image, '/', ...
        name_method_best_param, '/', ...
        'best_params/' ...   
    );
    
    % Get the contents of the directory
    contents = dir(save_best_folder_name);
    
    for idx_content = 1:length(contents)
        % Exclude '.' and '..' and check if the item is a folder
        if contents(idx_content).isdir && ~strcmp(contents(idx_content).name, '.') && ~strcmp(contents(idx_content).name, '..')
            name_params_savetext = append('best_params/' , contents(idx_content).name);
            break; % Exit the loop since there is only one folder
        end
    end

end


save_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_restored' ...
);


image_restored = HSI_restored(:,:,save_band)*restoration_Boost;

if ~isempty(arrow_methods_idc{1})
    for i = 1:numel(arrow_methods_idc)
        if find(arrow_methods_idc{i} == idx_method)
            image_restored = Embed_arrow(image_restored, ...
                        arrow_head_pos{i}, arrow_length, ...
                        arrow_handle_width, arrow_head_width, arrow_dir_tblr{i});
        end
    end
end

image_restored = Crop_Embed_image(image_restored, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

images_restored(:, :, idx_method) = image_restored;

end


%% Showing cropped images
if exist('is_plot', 'var') && is_plot
cat_restored_images = image_noisy;

for idx_method = idc_methods
    cat_restored_images = cat(2, cat_restored_images, images_restored(:,:,idx_method));
end

figure;
imshow(cat_restored_images);


end


%% Saving cropped images
if exist('is_save', 'var') && is_save
save_folder_name_cropped_image = append(...
    './result_image_for_JSTARS/', ...
    image ...
);

save_folder_name_restored_image = append(...
    save_folder_name_cropped_image, '/restored_image/', ...
     'b', num2str(save_band), '_m', num2str(restoration_Boost) ...
);


mkdir(save_folder_name_restored_image);


%% Saving each image
imwrite(image_noisy, append(save_folder_name_restored_image, '/image_noisy.png'));

for idx_method = idc_methods
    name_method = names_methods{idx_method};
    image_restored = images_restored(:,:,idx_method);
    
    imwrite(image_restored, ...
        append(save_folder_name_restored_image, ...
            '/image_', name_method, '.png') ...
    );

end
end
