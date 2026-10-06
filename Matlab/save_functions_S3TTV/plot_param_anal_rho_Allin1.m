clear;
close all;

addpath(genpath('sub_functions'))


%% Setting paramters
% boundary_width = 0;
% boundary_width = 2;
boundary_width = 3;
% boundary_width = 5;


is_save = 1;


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
% HSI_clean = single(HSI_clean(:, :, boundary_width + 1:end-boundary_width));
% HSI_noisy = single(HSI_noisy(:, :, boundary_width + 1:end-boundary_width));


name_method = 'S3TTV';


%% Setting common parameters
% boundary = 'Neumann';
boundary = 'circulant';

blocksize = [10,10]; % block size of STV-type methods
shiftstep = [1,1]; % [1,1] means full overlap

stopcri_index = 5;
stopcri = 10 ^ -stopcri_index;

maxiter = 20000;
% maxiter = 10;
% maxiter = 1;


%% Adjusting param
params_rho = {...
    0.91, ...
    0.92, ...
    0.93, ...
    0.94, ...
    0.95, ...
    0.96, ...
    0.97, ...
    0.98, ...
    0.99, ...
};

num_params_rho = numel(params_rho);


%% Preparing result list
mpsnr_list = zeros([num_params_rho, num_images]);
mssim_list = zeros([num_params_rho, num_images]);
% sam_list = zeros([num_params_rho, num_images]);


%% Loading MPSNRs and MSSIMs
for idx_image = idc_images
image = images{idx_image};

for idx_params_rho = 1:num_params_rho
rho = params_rho{idx_params_rho};

name_params_savetext = ...
        append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
            '_r', num2str(rho), '_stop1e-', num2str(stopcri_index)); % S3TTV

switch image
    case {'JasperRidge', 'PaviaU120'}
        save_folder_name = append(...
            './result_Param_Anal/Rho/' , ...
            'denoising_', image, '/', ...
            'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
                '_tint', num2str(deg.stripe_intensity), '/', ...
            name_method, '/', ...
            name_params_savetext, '/' ...   
        );
    case 'Beltsville'
        save_folder_name = append(...
            './result_Param_Anal/Rho/' , ...
            'denoising_', image, '/', ...
            'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
                '_pt', num2str(deg.stripe_rate), '/', ...
            name_method, '/', ...
            name_params_savetext, '/' ...   
        );
end


load(append(save_folder_name, 'image_result.mat'), ...
    'HSI_clean', 'HSI_restored' ...
);

HSI_restored = single(HSI_restored(:, :, boundary_width + 1:end-boundary_width));
HSI_clean = single(HSI_clean(:, :, boundary_width + 1:end-boundary_width));

mpsnr = calc_MPSNR(HSI_restored, HSI_clean);
mssim = calc_MSSIM(HSI_restored, HSI_clean);
% sam = calc_SAM(HSI_restored, HSI_clean);


mpsnr_list(idx_params_rho, idx_image) = mpsnr;
mssim_list(idx_params_rho, idx_image) = mssim;
% sam_list(idx_params_rho, idx_image) = sam;


end

end


%% Plotting and saving result
plot_style = {'LineWidth', 5.0};
label_style = {'FontSize', 100};
figure_style = {'Position', [100, 100, 1000, 800]};
ax_fontsize = 30;

index_rho = cell2mat(params_rho);

% figure();
figure(figure_style{:});
hold on
for idx_image = idc_images
    image_legend = images_legend{idx_image};
    plot(index_rho, mpsnr_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end
xlabel('\rho', label_style{:});
ylabel('MPSNR', label_style{:});
xlim([0.91, 0.99]);
ylim([28,37.5]);
ax = gca;
ax.FontSize = ax_fontsize;
xticks(index_rho);
xticklabels(arrayfun(@num2str, index_rho, 'UniformOutput', false));
legend('show')
grid on

if exist('is_save', 'var') && is_save == 1
    save_folder_name = append('./result_Param_Anal/Plot/');
    mkdir(save_folder_name)
    % saveas(gcf, append(save_file_name, 'rho_mpsnr.png'))
    exportgraphics(gcf, append(save_folder_name, 'rho_mpsnr.eps'), 'ContentType', 'image')
end


% figure();
figure(figure_style{:});
hold on
for idx_image = idc_images
    image_legend = images_legend{idx_image};
    plot(index_rho, mssim_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end
xlabel('\rho', label_style{:});
ylabel('MSSIM', label_style{:});
xlim([0.91, 0.99]);
ylim([0.76, 0.97]);
ax = gca;
ax.FontSize = ax_fontsize;
xticks(index_rho);
xticklabels(arrayfun(@num2str, index_rho, 'UniformOutput', false));
legend('show')
grid on
hold off

if exist('is_save', 'var') && is_save == 1
    % saveas(gcf, append(save_folder_name, 'rho_mssim.png'))
    exportgraphics(gcf, append(save_folder_name, 'rho_mssim.eps'), 'ContentType', 'image')
end


% figure(figure_style{:});
% plot(index_rho, sam_list, plot_style{:});
% xlabel('\rho', label_style{:});
% ylabel('sam', label_style{:});
% ax = gca;
% ax.FontSize = ax_fontsize;
% grid on
% saveas(gcf, append(save_file_name, 'rho_sam.png'))
