clear;
close all;

addpath(genpath('sub_functions'))
addpath('./func_metrics/')


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

names_methods_legend = { ...
    'SSTV', ...
    'HSSTV1', ...
    'HSSTV2', ...
    '$\ell_{0}$-$\ell_{1}$HTV', ...
    'STV', ...
    ...
    'SSST', ...
    'LRTDTV', ...
    'FGSLR', ...
    'TPTV', ...
    'FastHyMix', ...
    ...
    '$\mathrm{S_{3}TTV}$', ...
};


idc_methods = 1:numel(names_methods);
% idc_methods = [9, 10, 11];


idc_methods_best_param = 8:10;


is_show_spectra = 1;
% is_show_HSI = 1;

is_save_spectra = 1;


% is_show_sam_map = 1;

% idc_methods_for_sam_map = [1:numel(names_methods)];
idc_methods_for_sam_map = [2, numel(names_methods)];


% is_save_sam_map = 1;
% is_save_sam_map_color = 1;


%% Selecting conditions
noise_conditions = { ...
    {0.1,   0,      0,      0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,   0,      0},     ... % g0.05 ps0.05 pt0
    {0.1,   0.05,   0,      0},     ... % g0.1 ps0.05 pt0
    {0,     0,      0.05,   0.5},   ... % g0 ps0 pt0.05
    {0.05,  0,      0.05,   0.5},   ... % g0.05 ps0 pt0.05
    {0.1,   0,      0.05,   0.5},   ... % g0.1 ps0 pt0.05
    {0.05,  0.05,   0.05,   0.5},   ... % g0.05 ps0.05 pt0.05
    {0.1,   0.05,   0.05,   0.5},   ... % g0.1 ps0.05 pt0.05
};

idx_noise_condition = 8;

images = {...
    'JasperRidge', ...
    'PaviaU120', ...
    'Beltsville', ...
};

idx_image = 2;


switch idx_image
    case 1
        select_pos = [74, 92]; % Jasper Ridge

    case 2
        select_pos = [50, 89]; % Pavia Univerity tear1
        % select_pos = [74, 59]; % Pavia Univerity tear2
        % select_pos = [21, 103]; % Pavia Univerity tear2

    case 3
        % select_pos = [11, 84]; % Beltsville
        select_pos = [24, 17]; % Beltsville tear1
        % select_pos = [11, 75]; % Beltsville
end


%% Selecting noise condition
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
image = images{idx_image};


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


diff_magnification = 7;


%% Setting common parameters
if idx_noise_condition == 1
    rho = 0.98;
else
    rho = 0.95;
end


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
    ...
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)), ... % SSST
    append('stop1e-', num2str(stopcri_index)), ... % LRTDTV
    {}, ... % FGSLR
    {}, ... % TPTV
    {}, ... % FastHyMix
    ...
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

for idx_method = 1:length(contents)
    % Exclude '.' and '..' and check if the item is a folder
    if contents(idx_method).isdir && ~strcmp(contents(idx_method).name, '.') && ~strcmp(contents(idx_method).name, '..')
        names_params_savetext{idx_method_best_param} = append('best_params/' , contents(idx_method).name);
        break; % Exit the loop since there is only one folder
    end
end

end

%% Calculating clean and noisy HS image results
HSI_clean_spectrum = reshape(HSI_clean(select_pos(1), select_pos(2), :), [hsi.n3, 1]);


HSI_noisy_spectrum = reshape(HSI_noisy(select_pos(1), select_pos(2), :), [hsi.n3, 1]);
mpsnr  = calc_MPSNR(HSI_noisy, HSI_clean);
mssim  = calc_MSSIM(HSI_noisy, HSI_clean);
sam    = calc_SAM(HSI_noisy, HSI_clean);


fprintf('~~~ RESULTS ~~~\n');
fprintf('%s  \t MPSNR\t MSSIM\t SAM\n', blanks(max(strlength(names_methods))));
fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append('Noisy', blanks(max(strlength(names_methods)) - strlength('Noisy'))), ...
    mpsnr, mssim, sam);


%% Loading restored HS images
restored_spectra = zeros([numel(names_methods), hsi.n3]);
restored_sam_maps = zeros([hsi.n1, hsi.n2, numel(names_methods)]);
restored_HSIs = zeros([hsi.n1, hsi.n2, hsi.n3, numel(names_methods)]);


for idx_method = idc_methods
name_method = names_methods{idx_method};
name_params_savetext = names_params_savetext{idx_method};

load_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(load_folder_name, 'image_result.mat'), ...
    'HSI_restored' ...
);

load(append(load_folder_name, 'metric_vals.mat'), ...
    'mpsnr', 'mssim', 'sam' ...
);

[~, SAM_map] = calc_SAM(HSI_restored, HSI_clean);

fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
    append(name_method, blanks(max(strlength(names_methods)) - strlength(name_method))), ...
    mpsnr, mssim, sam);

% restored_spectra(idx_method, :) = reshape(HSI_restored(select_pos(1), select_pos(2), :), [1, hsi.n3]);
restored_spectra(idx_method, :) = HSI_restored(select_pos(1), select_pos(2), :);
restored_sam_maps(:, :, idx_method) = SAM_map;
restored_HSIs(:, :, :, idx_method) = HSI_restored;


end



%% Showing restored spectra
if exist('is_show_spectra', 'var') && is_show_spectra == 1
% colors = lines(numel(names_methods)); % 各手法用のカラーマップ
% colors = hsv(numel(names_methods)); % ユニークな色を生成するカラーマップ
colors = turbo(numel(names_methods));

figure('Position', [100, 100, 1500, 1500]); % [left, bottom, width, height]
hold on;

% Ground-truthのプロット
plot(HSI_clean_spectrum, ...
     'k', 'LineWidth', 2, 'DisplayName', 'Ground Truth'); % 赤色で太い線

% Noisy spectrumのプロット
% plot(HSI_noisy_spectrum, ...
%      'g--', 'LineWidth', 1.5, 'DisplayName', 'Noisy Spectrum'); % 緑の破線
plot(HSI_noisy_spectrum, ...
     'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Noisy Spectrum'); % 灰色の点線

% figure;
% plot(HSI_clean_spectrum, 'DisplayName', 'Ground truth', 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b')
% 
% hold on
% plot(HSI_noisy_spectrum, 'DisplayName', 'Noisy', 'LineWidth', 2, 'LineStyle', '--', 'Color', 'g')

for idx_method = idc_methods
    name_method = names_methods{idx_method};
    name_method_legend = names_methods_legend{idx_method};
    if idx_method == numel(names_methods)
        plot(restored_spectra(idx_method, :), ...
            'r', 'DisplayName', name_method_legend, 'LineWidth', 2)
    else
        plot(restored_spectra(idx_method, :), ...
            'Color', colors(idx_method, :), ...
             'DisplayName', name_method_legend, 'LineWidth', 2, 'LineStyle', '--');
    end
end

set(gca, 'FontSize', 35);

xlim([1, hsi.n3])
% ylim([min(HSI_clean_spectrum, [], 'all')*0.7, max(HSI_clean_spectrum, [], 'all')*1.2])
ylim([0, max(HSI_clean_spectrum, [], 'all')*1.5])
% プロットの設定
% legend('show', 'Location', 'best', 'Interpreter', 'latex', 'NumColumns', 2); % 凡例の表示
legend('show', 'Location', 'northwest', 'Interpreter', 'latex', 'NumColumns', 2); %Pavia University  
% legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'NumColumns', 2); % Beltsville
xlabel('Spectral Band'); % x軸ラベル
ylabel('Intensity'); % y軸ラベル
grid on;
hold off;

if exist('is_save_spectra', 'var') && is_save_spectra
    save_folder_name = append(...
        './result_image_for_JSTARS/' , ...
        'spectra/' ...
    );
    mkdir(save_folder_name);
    
    save_spectra_file_name = append(...
        save_folder_name, ...
        image, '_', ...
        'spectra_', ...
        num2str(select_pos(1)), '_', num2str(select_pos(2)), '.png' ...
    );

    saveas(gcf, save_spectra_file_name); 
end

end


%% Showing HSI
if exist('is_show_HSI', 'var') && is_show_HSI == 1
cat_HSI = cat(2, HSI_clean, HSI_noisy);
for idx_method = idc_methods
    cat_HSI = cat(2, cat_HSI, restored_HSIs(:,:,:,idx_method));
end

HSI_cleans = repmat(HSI_clean, [1, numel(idc_methods)+2, 1]);

cat_diff = abs(cat_HSI - HSI_cleans) * diff_magnification;
show_HSI = cat(1, cat_HSI, cat_diff);
implay(show_HSI)
end


%% Showing SAM maps
if exist('is_show_sam_map', 'var') && is_show_sam_map == 1
cat_SAM_map = [];

for idx_method = idc_methods_for_sam_map
    cat_SAM_map = cat(2, cat_SAM_map, restored_sam_maps(:,:,idx_method));
end

figure;
imagesc(cat_SAM_map)

end


%% Saving SAM maps
if exist('is_save_sam_map', 'var') && is_save_sam_map == 1
save_sam_map_folder_name = append(...
    './result_image_for_JSTARS/' , ...
    'SAM_map/', ...
    image, '/' ...
);

mkdir(save_sam_map_folder_name)


for idx_method = idc_methods
    name_method = names_methods{idx_method};
    sam_map = restored_sam_maps(:,:,idx_method);

    imwrite(sam_map, ...
        append(save_sam_map_folder_name, ...
        'sam_map_', name_method, '.png' ...
        ) ...
    );
end

end

%% Saving SAM maps colored
if exist('is_save_sam_map_color', 'var') && is_save_sam_map_color == 1
save_sam_map_folder_name = append(...
    './result_image_for_JSTARS/' , ...
    'SAM_map_color/', ...
    image, '/' ...
);

mkdir(save_sam_map_folder_name)

types = 'hot';
% types = 'gray';
% types = 'turbo';
% types = 'jet';
% types = 'parula';

cmap = colormap(types);
close all;


for idx_method = idc_methods
    name_method = names_methods{idx_method};
    sam_map_tmp = restored_sam_maps(:,:,idx_method);

    sam_map_gray = sam_map_tmp ./ 1.3;
    sam_map_color = Gray2RGB(sam_map_gray, cmap);
    image_restored = Crop_Embed_image(image_restored, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

    imwrite(sam_map_color, ...
        append(save_sam_map_folder_name, ...
        'sam_map_', name_method, '.png' ...
        ) ...
    );
end

end