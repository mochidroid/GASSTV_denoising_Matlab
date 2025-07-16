clear;
close all;

addpath(genpath('sub_functions'))
addpath('func_metrics')


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
% idc_methods = [2, 3, numel(names_methods)];


idc_methods_best_param = 8:10;


% boundary_width = 0;
% boundary_width = 2;
boundary_width = 3;
% boundary_width = 5;


% is_plot_metrics = 1;
% is_plot_LaTeX_table_mpsnr = 1;
% is_plot_LaTeX_table_mssim = 1;
is_plot_LaTeX_table_sam = 1;


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

num_cases = numel(noise_conditions);

idc_noise_conditions = 1:num_cases;
% idc_noise_conditions = [1, 3];
% idc_noise_conditions = 8;

images = {...
    'JasperRidge', ...
    'PaviaU120', ...
    'Beltsville', ...
};

images_outputs = {...
    'Jasper Ridge', ...
    'Pavia University', ...
    'Beltsville', ...
};

idc_images = 1:numel(images);
% idc_images = 1;



%% Loading HS images and Calculating metrics
vals_metrics = zeros([numel(names_methods), 3, numel(noise_conditions), numel(images)]);
vals_metrics_noisy = zeros([1, 3, numel(noise_conditions), numel(images)]);

for idx_image = idc_images
for idx_noise_condition = idc_noise_conditions
%% Selecting noise condition
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
image = images{idx_image};
image_output = images_outputs{idx_image};


[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean(:, :, boundary_width + 1:end-boundary_width));
HSI_noisy = single(HSI_noisy(:, :, boundary_width + 1:end-boundary_width));


%% Setting common parameters
switch idx_noise_condition
    case 1
        rho = 0.98;
    case 4
        rho = 0.98;
        epsilon_rho = 0.01;
    otherwise
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
save_S3TTV_Calc_FV_folder_name = append(...
    './result/' , ...
    'denoising_', image, '/', ...
    'g', num2str(deg.gaussian_sigma), '_ps', num2str(deg.sparse_rate), ...
        '_pt', num2str(deg.stripe_rate), '/', ...
    'S3TTV/Calc_FV' ...
);

switch idx_noise_condition
    case 4
        % Check if Calc_FV folder exists
        if isfolder(save_S3TTV_Calc_FV_folder_name)
            name_param_save_text_S3TTV = ...
                append('Calc_FV/', ...
                    'bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
                    '_r', num2str(rho), '_er', num2str(epsilon_rho), '_stop1e-', num2str(stopcri_index) ...
                );
        else
            name_param_save_text_S3TTV = ...
                append(...
                    'bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
                    '_r', num2str(rho), '_er', num2str(epsilon_rho), '_stop1e-', num2str(stopcri_index) ...
                );
        end

        names_params_savetext = { ...
            append('r', num2str(rho), '_er', num2str(epsilon_rho), '_stop1e-', num2str(stopcri_index)), ... % SSTV
            append('o', num2str(HSSTV_omega), '_r', num2str(rho), '_er', num2str(epsilon_rho), ...
                '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L1
            append('o', num2str(HSSTV_omega), '_r', num2str(rho), '_er', num2str(epsilon_rho), ...
                '_stop1e-', num2str(stopcri_index)), ... % HSSTV_L12
            append('sr', num2str(l0l1HTV_stepsize_reduction), ...
                '_r', num2str(rho), '_er', num2str(epsilon_rho), '_maxiter', num2str(maxiter)) ... % l0l1HTV
            append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
                '_r', num2str(rho), '_er', num2str(epsilon_rho), '_stop1e-', num2str(stopcri_index)), ... % STV
            append('bl', num2str(blocksize(1)), '_st', num2str(shiftstep(1)), ...
                '_r', num2str(rho), '_er', num2str(epsilon_rho), '_stop1e-', num2str(stopcri_index)), ... % SSST
            append('stop1e-', num2str(stopcri_index)), ... % LRTDTV
            {}, ... % FGSLR
            {}, ... % TPTV
            {}, ... % FastHyMix
            ...
            name_param_save_text_S3TTV, ... % S3TTV
        };

    otherwise
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
end

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

%% Loading restored HS images
mpsnr_noisy  = calc_MPSNR(HSI_noisy, HSI_clean);
mssim_noisy  = calc_MSSIM(HSI_noisy, HSI_clean);
sam_noisy    = calc_SAM(HSI_noisy, HSI_clean);

vals_metrics_noisy(1, :, idx_noise_condition, idx_image) = [mpsnr_noisy, mssim_noisy, sam_noisy];


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

HSI_restored = single(HSI_restored(:, :, boundary_width + 1:end-boundary_width));


% list_psnr_per_band = cat(1, list_psnr_per_band, psnr_per_band);
% list_ssim_per_band = cat(1, list_ssim_per_band, ssim_per_band);


mpsnr = calc_MPSNR(HSI_restored, HSI_clean);
mssim = calc_MSSIM(HSI_restored, HSI_clean);
sam = calc_SAM(HSI_restored, HSI_clean);

sam = sam * pi / 180;

vals_metrics(idx_method, :, idx_noise_condition, idx_image) = [mpsnr, mssim, sam];


end

end
end


%% Plotting metrics
if exist('is_plot_metrics', 'var') && is_plot_metrics == 1

for idx_image = idc_images
for idx_noise_condition = idc_noise_conditions
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    image_output = images_outputs{idx_image};

    fprintf('~~~ SETTINGS ~~~\n');
    fprintf('Image: %s Size: (%d, %d, %d)\n', image_output, hsi.n1, hsi.n2, hsi.n3);
    fprintf('Case %d: g%g ps%g pt%g\n\n', idx_noise_condition, deg.gaussian_sigma, deg.sparse_rate, deg.stripe_intensity);

    fprintf('~~~ RESULTS ~~~\n');
    fprintf('%s  \t MPSNR\t MSSIM\t SAM\n', blanks(max(strlength(names_methods))));
    
    fprintf('%s: \t %#.4g\t %#.4g\t %#.4g\n', ...
        append('noisy', blanks(max(strlength(names_methods)) - 5)), ...
        mpsnr_noisy, mssim_noisy, sam_noisy);
    
    
    for idx_method = idc_methods
        name_method = names_methods{idx_method};
    
        fprintf('%s: \t %g\t %g\t %g\n', ...
            append(name_method, blanks(max(strlength(names_methods)) - strlength(name_method))), ...
            vals_metrics(idx_method, 1), vals_metrics(idx_method, 2), vals_metrics(idx_method, 3));
    end

end
end

end


%% Plotting LaTeX tables for MPSNR
if exist('is_plot_LaTeX_table_mpsnr', 'var') && is_plot_LaTeX_table_mpsnr == 1
name_folder_o_tex = './result_image_for_JSTARS/plot_LaTeX_table';
mkdir(name_folder_o_tex);
name_file_o = append(name_folder_o_tex, '/table_mpsnr_bo', num2str(boundary_width), '.txt');
fid = fopen(name_file_o, 'w');


text_tex = '';

for idx_image = idc_images
    image_output = images_outputs{idx_image};

    text_tex = append(text_tex, '\\cmidrule(lr){1-13} \n\n');
    text_tex = append(text_tex, '\\multirow{', num2str(num_cases), '}{*}{', image_output, '} \n');

    for idx_noise_condition = idc_noise_conditions
        text_tex = append(text_tex, '& Case ', num2str(idx_noise_condition), ' & \n');

        % Stack first and second best
        vals_mpsnr_tmp = vals_metrics(:, 1, idx_noise_condition, idx_image);
        
        [vals_mpsnr_tmp_sorted, idx_mpsnr_tmp_sorted] = sort(vals_mpsnr_tmp, "descend");
        idx_first_second_best = zeros(size(vals_mpsnr_tmp));
        num_best = sum(double(vals_mpsnr_tmp_sorted == vals_mpsnr_tmp_sorted(1)), "all");
        idx_first_second_best(idx_mpsnr_tmp_sorted(1:num_best)) = 1;
        num_second = sum(double(vals_mpsnr_tmp_sorted == vals_mpsnr_tmp_sorted(num_best + 1)), "all");
        idx_first_second_best(idx_mpsnr_tmp_sorted((1 + num_best):(num_second + num_best))) = 2;

        for idx_method = idc_methods
            mpsnr_text_tmp = sprintf('%#.4g', vals_metrics(idx_method, 1, idx_noise_condition, idx_image));

            if idx_first_second_best(idx_method) == 1
                mpsnr_text_tmp = append('\\textbf{', mpsnr_text_tmp, '}');
            elseif idx_first_second_best(idx_method) == 2
                mpsnr_text_tmp = append('\\underline{', mpsnr_text_tmp, '}');
            end
            text_tex = append(text_tex, mpsnr_text_tmp);
            if idx_method ~= idc_methods(end)
                text_tex = append(text_tex, ' & ');
            elseif idx_method == idc_methods(end)
                text_tex = append(text_tex, ' \\\\ \n');
            end
        end
    
    end

    text_tex = append(text_tex, '\n');

    if idx_image == idc_images(end)
        text_tex = append(text_tex, '\\bottomrule');
    end
end

fprintf(fid, text_tex);
end


%% Plotting LaTeX tables for MSSIM
if exist('is_plot_LaTeX_table_mssim', 'var') && is_plot_LaTeX_table_mssim == 1
name_folder_o_tex = './result_image_for_JSTARS/plot_LaTeX_table';
mkdir(name_folder_o_tex);
name_file_o = append(name_folder_o_tex, '/table_mssim_bo', num2str(boundary_width), '.txt');
fid = fopen(name_file_o, 'w');


text_tex = '';

for idx_image = idc_images
    image_output = images_outputs{idx_image};

    text_tex = append(text_tex, '\\cmidrule(lr){1-13} \n\n');
    text_tex = append(text_tex, '\\multirow{', num2str(num_cases), '}{*}{', image_output, '} \n');

    for idx_noise_condition = idc_noise_conditions
        text_tex = append(text_tex, '& Case ', num2str(idx_noise_condition), ' & \n');

        % Stack first and second best
        vals_mssim_tmp = vals_metrics(:, 2, idx_noise_condition, idx_image);
        [vals_mssim_tmp_sorted, idx_mssim_tmp_sorted] = sort(vals_mssim_tmp, "descend");
        idx_first_second_best = zeros(size(vals_mssim_tmp));
        num_best = sum(double(vals_mssim_tmp_sorted == vals_mssim_tmp_sorted(1)), "all");
        idx_first_second_best(idx_mssim_tmp_sorted(1:num_best)) = 1;
        num_second = sum(double(vals_mssim_tmp_sorted == vals_mssim_tmp_sorted(num_best + 1)), "all");
        idx_first_second_best(idx_mssim_tmp_sorted((1 + num_best):(num_second + num_best))) = 2;

        for idx_method = idc_methods
            mssim_text_tmp = sprintf('%#.4g', vals_metrics(idx_method, 2, idx_noise_condition, idx_image));

            if idx_first_second_best(idx_method) == 1
                mssim_text_tmp = append('\\textbf{', mssim_text_tmp, '}');
            elseif idx_first_second_best(idx_method) == 2
                mssim_text_tmp = append('\\underline{', mssim_text_tmp, '}');
            end
            text_tex = append(text_tex, mssim_text_tmp);
            if idx_method ~= idc_methods(end)
                text_tex = append(text_tex, ' & ');
            elseif idx_method == idc_methods(end)
                text_tex = append(text_tex, ' \\\\ \n');
            end
        end
    
    end

    text_tex = append(text_tex, '\n');

    if idx_image == idc_images(end)
        text_tex = append(text_tex, '\\bottomrule');
    end

end

fprintf(fid, text_tex);
end


%% Plotting LaTeX tables for SAM
if exist('is_plot_LaTeX_table_sam', 'var') && is_plot_LaTeX_table_sam == 1
name_folder_o_tex = './result_image_for_JSTARS/plot_LaTeX_table';
mkdir(name_folder_o_tex);
% name_file_o = append(name_folder_o_tex, '/table_sam_bo', num2str(boundary_width), '.txt');
name_file_o = append(name_folder_o_tex, '/table_sam_bo', num2str(boundary_width), '_degree.txt');
fid = fopen(name_file_o, 'w');


text_tex = '';

for idx_image = idc_images
    image_output = images_outputs{idx_image};

    text_tex = append(text_tex, '\\cmidrule(lr){1-13} \n\n');
    text_tex = append(text_tex, '\\multirow{', num2str(num_cases), '}{*}{', image_output, '} \n');

    for idx_noise_condition = idc_noise_conditions
        text_tex = append(text_tex, '& Case ', num2str(idx_noise_condition), ' & \n');

        % Stack first and second best
        vals_sam_tmp = vals_metrics(:, 3, idx_noise_condition, idx_image);
        [vals_sam_tmp_sorted, idx_sam_tmp_sorted] = sort(vals_sam_tmp, "ascend");
        idx_first_second_best = zeros(size(vals_sam_tmp));
        num_best = sum(double(vals_sam_tmp_sorted == vals_sam_tmp_sorted(1)), "all");
        idx_first_second_best(idx_sam_tmp_sorted(1:num_best)) = 1;
        num_second = sum(double(vals_sam_tmp_sorted == vals_sam_tmp_sorted(num_best + 1)), "all");
        idx_first_second_best(idx_sam_tmp_sorted((1 + num_best):(num_second + num_best))) = 2;

        for idx_method = idc_methods
            sam_text_tmp = sprintf('%#.4g', vals_metrics(idx_method, 3, idx_noise_condition, idx_image));

            if idx_first_second_best(idx_method) == 1
                sam_text_tmp = append('\\textbf{', sam_text_tmp, '}');
            elseif idx_first_second_best(idx_method) == 2
                sam_text_tmp = append('\\underline{', sam_text_tmp, '}');
            end
            text_tex = append(text_tex, sam_text_tmp);
            if idx_method ~= idc_methods(end)
                text_tex = append(text_tex, ' & ');
            elseif idx_method == idc_methods(end)
                text_tex = append(text_tex, ' \\\\ \n');
            end
        end
    
    end

    text_tex = append(text_tex, '\n');

    if idx_image == idc_images(end)
        text_tex = append(text_tex, '\\bottomrule');
    end

end

fprintf(fid, text_tex);
end