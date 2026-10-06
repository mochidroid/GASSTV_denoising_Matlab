clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Switching operator
is_plot_metrics = 1;
is_show_HSI = 1;
diff_magnification = 7;

is_plot_psnr_and_ssim_per_band = 1;
% is_output_csv = 1;

num_show = 10;


%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % g0.05 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.05,  0,     0.05,  0.5,  0},     ... % g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
    {0.1,   0,     0,     0,    0.01},  ... % g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % g0.1 ps0.05 pt0.05 pd0.01
};
% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = 7;

images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

% images = {...
%     "JasperRidge64", ...
%     "PaviaU64", ...
% };

% idc_images = 1:numel(images);
idc_images = [1:2];
% idc_images = 3;


%% Setting common parameters
rhos = {0.95};

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;
% maxiter = 5;

load("dir_save_comp_folder.mat", "dir_save_comp_folder");


%% Setting each methods info
% GASSTV_CondatVu_OG
% params_GASSTV.sigma_sp = [0.01, 0.1, 1];
params_GASSTV.sigma_sp = {"90", "med"};
params_GASSTV.sigma_l = {"90", "med"};
params_GASSTV.lambda_rho_sp = [0.9, 1, 1.1];
% params_GASSTV.sigma_s = [0.01, 1];
% params_GASSTV.lambda1 = [0.01, 0.1, 1];
% params_GASSTV.lambda1 = [0.01];
params_GASSTV.lambda2 = [0.01, 0.1, 1, 10];
params_GASSTV.num_segments = [3, 5, 10];
% params_GASSTV.num_segments = [5];
% params_GASSTV.order_filt = [1, 3];
params_GASSTV.order_filt = [1];
methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu_OG", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{params_GASSTV.lambda_rho_sp, params_GASSTV.lambda2, params_GASSTV.sigma_sp, params_GASSTV.sigma_l, ...
        params_GASSTV.num_segments, params_GASSTV.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", true ...
);


% GASSTV_CondatVu
% params_GASSTV.sigma_sp = [0.01, 0.1, 1];
params_GASSTV.sigma_sp = {"90", "med"};
params_GASSTV.sigma_l = {"90", "med"};
params_GASSTV.lambda_rho_sp = [0.9, 1, 1.1];
params_GASSTV.lambda2 = [0.1, 1, 10];
params_GASSTV.num_segments = [3, 5, 10];
% params_GASSTV.num_segments = [5];
params_GASSTV.order_filt = [5, 1];
% params_GASSTV.order_filt = [5];
methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{params_GASSTV.lambda_rho_sp, params_GASSTV.lambda2, params_GASSTV.sigma_sp, params_GASSTV.sigma_l, ...
        params_GASSTV.num_segments, params_GASSTV.order_filt, maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", true ...
);


methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Loading each method
for idx_method = 1:num_methods
    name_method = methods_info(idx_method).name;
    params_name = methods_info(idx_method).param_names;
    params_cell = methods_info(idx_method).params;
    
    % Generate all parameter combinations
    [params_comb, num_params_comb] = ParamsList2Comb(params_cell);

    %% Loading each condition
    for idx_noise_condition = idc_noise_conditions
        for idx_image = idc_images
            
            %% Generating observation info
            deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
            deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
            deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
            deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
            deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
            image = images{idx_image};
            
            dir_result_folder = fullfile(...
                dir_save_comp_folder, ...
                append("denoising_", image), ...
                append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                    "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
                name_method ...
            );
            
            fprintf("\n================================================================\n");
            fprintf("Analyzing Method: %s\n", name_method);
            fprintf("Image: %s (g%.2f, ps%.2f, pt%.2f)\n", image, deg.gaussian_sigma, deg.sparse_rate, deg.stripe_rate);
            fprintf("Total Parameter Combinations: %d\n", num_params_comb);
            fprintf("================================================================\n");

            %% 1. Load All Results into Struct Array
            results_list = struct();
            valid_count = 0;
            
            for idx_params_comb = 1:num_params_comb
                % Construct params struct from combination
                params = struct();
                for idx_params = 1:numel(params_name)
                    params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
                end
                
                % Generate filename using the defined function
                name_params_savetext = methods_info(idx_method).get_params_savetext(params);
                file_path = fullfile(dir_result_folder, append(name_params_savetext, ".mat"));
                
                % Try loading
                if isfile(file_path)
                    try
                        data = load(file_path, "val_mpsnr", "val_mssim", "val_sam");
                        valid_count = valid_count + 1;
                        
                        % Metrics
                        results_list(valid_count).MPSNR = data.val_mpsnr;
                        results_list(valid_count).MSSIM = data.val_mssim;
                        results_list(valid_count).SAM   = data.val_sam;
                        
                        % Parameters (Explicitly store for analysis)
                        results_list(valid_count).lambda_rho_sp = params.lambda_rho_sp;
                        results_list(valid_count).lambda2       = params.lambda2;
                        results_list(valid_count).sigma_sp      = string(params.sigma_sp);
                        results_list(valid_count).sigma_l       = string(params.sigma_l);
                        results_list(valid_count).num_segments  = params.num_segments;
                        results_list(valid_count).order_filt    = params.order_filt;
                        results_list(valid_count).NameText      = string(name_params_savetext);
                        
                    catch
                        fprintf("Error loading: %s\n", name_params_savetext);
                    end
                else
                    % fprintf("File not found: %s\n", name_params_savetext);
                end
            end
            
            if valid_count == 0
                fprintf("No valid results loaded for this condition.\n");
                continue;
            end
            
            % Convert to Table
            T = struct2table(results_list);
            
            %% 2. Display Best Results
            fprintf("\n--- Top %d Best MPSNR Settings ---\n", num_show);
            T_sorted = sortrows(T, 'MPSNR', 'descend');
            disp(T_sorted(1:min(num_show, height(T)), {'MPSNR', 'MSSIM', 'lambda_rho_sp', 'lambda2', 'sigma_sp', 'sigma_l', 'num_segments'}));
            
            %% 3. Visualization
            
            % Get Best Parameters to fix others while plotting heatmaps
            best_row = T_sorted(1, :);
            best_sig_sp = best_row.sigma_sp;
            best_sig_l  = best_row.sigma_l;
            best_seg    = best_row.num_segments;
            best_filt   = best_row.order_filt;
            
            % --- (A) Heatmap: Lambda2 vs Lambda_rho_sp ---
            % Filter data matching best discrete params
            mask_hm = (T.sigma_sp == best_sig_sp) & ...
                      (T.sigma_l == best_sig_l) & ...
                      (T.num_segments == best_seg) & ...
                      (T.order_filt == best_filt);
            T_hm = T(mask_hm, :);
            
            if height(T_hm) > 1
                figure('Name', append(image, ": Lambda Heatmap"));
                
                x_vals = unique(T_hm.lambda2);
                y_vals = unique(T_hm.lambda_rho_sp);
                Z_mat = nan(length(y_vals), length(x_vals));
                
                for i = 1:length(y_vals)
                    for j = 1:length(x_vals)
                        row = T_hm(T_hm.lambda_rho_sp == y_vals(i) & T_hm.lambda2 == x_vals(j), :);
                        if ~isempty(row)
                            Z_mat(i, j) = row.MPSNR(1);
                        end
                    end
                end
                
                imagesc(Z_mat);
                colorbar;
                title({append("MPSNR Heatmap: ", image), ...
                       sprintf("Fixed: sig_{sp}=%s, sig_{l}=%s, seg=%d", best_sig_sp, best_sig_l, best_seg)}, ...
                       'Interpreter', 'none');
                xlabel('\lambda_2 (Spectral)');
                ylabel('\lambda_{\rho, sp} (Spatial)');
                set(gca, 'XTick', 1:length(x_vals), 'XTickLabel', string(x_vals));
                set(gca, 'YTick', 1:length(y_vals), 'YTickLabel', string(y_vals));
                axis xy;
            end
            
            % --- (B) Boxplots: Discrete Parameters Impact ---
            figure('Name', append(image, ": Parameter Impacts"), 'Position', [100, 100, 1400, 400]);
            
            % Subplot 1: Sigma Combinations
            subplot(1, 4, 1);
            group_sig = strcat("sp:", T.sigma_sp, newline, "l:", T.sigma_l);
            boxplot(T.MPSNR, group_sig);
            title('Sigma Settings');
            ylabel('MPSNR');
            grid on;
            
            % Subplot 2: Num Segments
            subplot(1, 4, 2);
            boxplot(T.MPSNR, T.num_segments);
            title('Num Segments');
            xlabel('K');
            grid on;
            
            % Subplot 3: Filter Order
            subplot(1, 4, 3);
            boxplot(T.MPSNR, T.order_filt);
            title('Filter Order');
            xlabel('Order');
            grid on;

            % Subplot 4: Lambda2 (Overall Distribution)
            subplot(1, 4, 4);
            boxplot(T.MPSNR, T.lambda2);
            title('Lambda2 Overall');
            xlabel('\lambda_2');
            grid on;
            
            drawnow;
        end
    end
end
fprintf("******* Analysis Finished *******\n");