clear;
close all;

addpath('./func_metrics')
addpath(genpath('sub_functions'))


%% Setting paramters
% boundary_width = 0;
% boundary_width = 2;
boundary_width = 3;
% boundary_width = 5;


% is_save = 1;


%% Generating observation
deg.gaussian_sigma      = 0.1;
deg.sparse_rate         = 0.05;
deg.stripe_rate         = 0.05;
deg.stripe_intensity    = 0.5;

images = {...
    'JasperRidge', ...
    'PaviaU120', ...
    'Beltsville', ...
};


images_legend = { ...
    'Jasper Ridge', ...
    'Pavia University', ...
    'Beltsville', ...
};

idc_images = 1:numel(images);
% idc_images = 1;

num_images = numel(idc_images);


% [HSI_clean, hsi] = Load_HSI(image);
% [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);
% 
% HSI_clean = single(HSI_clean);
% HSI_noisy = single(HSI_noisy);


name_method = 'S3TTV';


%% Setting common parameters
% boundary = 'Neumann';
boundary = 'circulant';

% blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

rho = 0.95;

stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;

maxiter = 20000;
% maxiter = 10;
% maxiter = 1;


%% Adjusting param
params_blocksize = {...
    [2, 2], ...
    [4, 4], ...
    [6, 6], ...
    [8, 8], ...
    [10, 10], ...
    [12, 12], ...
    [14, 14], ...
    [16, 16], ...
};

num_params_blocksize = numel(params_blocksize);


mpsnr_list = zeros([num_params_blocksize, num_images]);
mssim_list = zeros([num_params_blocksize, num_images]);
% sam_list = zeros([num_params_blocksize, num_images]);

%% Loading MPSNRs and MSSIMs
for idx_image = idc_images
image = images{idx_image};

for idx_params_blocksize = 1:num_params_blocksize
blocksize = params_blocksize{idx_params_blocksize};

name_params_savetext = ...
    append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
        '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)); % S3TTV

save_folder_name = append(...
    './result_Param_Anal/Blocksize/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_tint', num2str(deg.stripe_intensity), '/', ...
    name_method, '/', ...
    name_params_savetext, '/' ...   
);

load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_clean', 'HSI_restored' ...
);

HSI_restored = single(HSI_restored(:, :, boundary_width + 1:end-boundary_width));
HSI_clean = single(HSI_clean(:, :, boundary_width + 1:end-boundary_width));

mpsnr = calc_MPSNR(HSI_restored, HSI_clean);
mssim = calc_MSSIM(HSI_restored, HSI_clean);
% sam = calc_SAM(HSI_restored, HSI_clean);


mpsnr_list(idx_params_blocksize, idx_image) = mpsnr;
mssim_list(idx_params_blocksize, idx_image) = mssim;
% sam_list(idx_params_blocksize, idx_image) = sam;


end

end

%% Plotting and saving result
% plot_style = {'LineWidth', 2.0};
plot_style = {'LineWidth', 5.0};
label_style = {'FontSize', 25};
% label_style = {'FontSize', 300};
% figure_style = {'Position', [100, 100, 800, 600]};
figure_style = {'Position', [100, 100, 600, 500]};
ax_fontsize = 17;
blocksize_labels = cellfun(@(x) sprintf('[%d,%d]', x(1), x(2)), params_blocksize, 'UniformOutput', false);

index_blocksize = cell2mat(params_blocksize);

figure(figure_style{:});
hold on
for idx_image = idc_images
    image_legend = images_legend{idx_image};
    plot(1:num_params_blocksize, mpsnr_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end


xlabel('Block size', label_style{:});
ylabel('MPSNR', label_style{:});
xlim([1, 8]);
ylim([33, 37]);
ax = gca;
ax.FontSize = ax_fontsize;
xticks(1:length(blocksize_labels));
xticklabels(blocksize_labels);
grid on
legend('show')

hold off

if exist('is_save', 'var') && is_save == 1
    save_folder_name = append('./result_Param_Anal/Plot/');
    mkdir(save_folder_name)
    % saveas(gcf, append(save_folder_name, 'blocksize_mpsnr.png'))
    exportgraphics(gcf, append(save_folder_name, 'blocksize_mpsnr.eps'), 'ContentType', 'image')
end


figure(figure_style{:});
hold on
for idx_image = idc_images
    disp(idx_image)
    image_legend = images_legend{idx_image};
    plot(1:num_params_blocksize, mssim_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end
xlabel('Block size', label_style{:});
ylabel('MSSIM', label_style{:});
xlim([1, 8]);
% ylim([0.91, 0.93]) %JasperRidge
% ylim([0.87, 0.90]) %PaviaU120
ylim([0.87, 0.94]);
ax = gca;
ax.FontSize = ax_fontsize;
xticks(1:length(blocksize_labels));
xticklabels(blocksize_labels);
grid on
legend('show')

if exist('is_save', 'var') && is_save == 1
    % saveas(gcf, append(save_file_name, 'blocksize_mssim.png'))
    exportgraphics(gcf, append(save_folder_name, 'blocksize_mssim.eps'), 'ContentType', 'image')

end


% figure(figure_style{:});
% plot(1:length(sam_list), sam_list, plot_style{:});
% xlabel('Block size', label_style{:});
% ylabel('sam', label_style{:});
% ax = gca;
% ax.FontSize = ax_fontsize;
% xticks(1:length(blocksize_labels));
% xticklabels(blocksize_labels);
% grid on
