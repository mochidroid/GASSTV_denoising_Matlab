clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

fprintf("******* Simulation Data Parameter Search (MPSNR/MSSIM/SAM) *******\n");

%% 1. Config & Data Selection
% =========================================================================
% 表示する上位件数
num_show = 10; 

% ノイズ条件の設定 (シミュレーション用)
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % 1. Gaussian only
    {0.05,  0.05,  0,     0,    0},     ... % 2. Mixed 1
    {0.1,   0.05,  0,     0,    0},     ... % 3. Mixed 2
    {0.05,  0,     0.05,  0.5,  0},     ... % 4. Stripe 1
    {0.1,   0,     0.05,  0.5,  0},     ... % 5. Stripe 2
    {0.05,  0.05,  0.05,  0.5,  0},     ... % 6. Complex 1
    {0.1,   0.05,  0.05,  0.5,  0},     ... % 7. Complex 2 (Target)
};
% 評価したいノイズ条件のインデックス
idc_noise_conditions = 7; 

% 画像の設定
images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};
idc_images = 1; % 1:JasperRidge

% 結果保存フォルダのロード
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

% 共通パラメータ
maxiter = 20000;
stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;
rhos = {0.95};
fmt4s = @(x) round(x,4,"significant");


%% 2. Setting Method Parameters
% =========================================================================
% GASSTV_CondatVu パラメータ設定
GASSTV_CondatVu.sigma_sp = {"90", "med"};
GASSTV_CondatVu.sigma_l = {"med"};
GASSTV_CondatVu.lambda_rho_sp = [0.9, 1, 1.1];
GASSTV_CondatVu.lambda2 = [0.1, 1, 10];
GASSTV_CondatVu.num_segments = [4];
GASSTV_CondatVu.order_filt = [5, 1];

methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", "rho_radius"}}, ...
    "params", {{GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, ...
        GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, ...
        maxiter, stopcri, rhos}}, ...
    "get_params_savetext", @(params) ...
        sprintf("l%.2g_%.2g_sig%s_%s_ns%d_fil%d_r%.2f_stop1e-%d", ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, params.rho_radius, stopcri_idx), ...
    "enable", true ...
);

methods_info = methods_info([methods_info.enable]); 
num_methods = numel(methods_info);


%% 3. Main Loop: Load Metrics -> Sort -> Save Best
% =========================================================================
for idx_method = 1:num_methods
    name_method = methods_info(idx_method).name;
    params_name = methods_info(idx_method).param_names;
    params_cell = methods_info(idx_method).params;
    
    % パラメータ組み合わせの生成
    [params_comb, num_params_comb] = ParamsList2Comb(params_cell);
    
    for idx_noise = idc_noise_conditions
        for idx_img = idc_images
            
            % --- 条件設定 ---
            deg.gaussian_sigma   = noise_conditions{idx_noise}{1};
            deg.sparse_rate      = noise_conditions{idx_noise}{2};
            deg.stripe_rate      = noise_conditions{idx_noise}{3};
            deg.deadline_rate    = noise_conditions{idx_noise}{5};
            image_name = images{idx_img};
            
            % 結果フォルダパス
            dir_result_folder = fullfile(...
                dir_save_comp_folder, ...
                append("denoising_", image_name), ...
                append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                    "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
                name_method ...
            );
            
            fprintf('\n~~~ %s | %s | Case %d ~~~\n', name_method, image_name, idx_noise);
            
            % --- 結果収集ループ ---
            summary_metrics = struct([]);
            
            for idx_comb = 1:num_params_comb
                % パラメータ構造体の作成
                params = struct();
                for k = 1:numel(params_name)
                    params.(params_name{k}) = params_comb{idx_comb}{k};
                end
                
                % ファイル名取得
                param_savetext = methods_info(idx_method).get_params_savetext(params);
                
                mat_file = fullfile(dir_result_folder, append(param_savetext, ".mat"));
                
                if isfile(mat_file)
                    % 結果ロード (val_mpsnr, val_mssim, val_sam)
                    load(mat_file, "val_mpsnr", "val_mssim", "val_sam");
                    
                    % 構造体に格納
                    metrics_struct = params;
                    metrics_struct.MPSNR = fmt4s(val_mpsnr);
                    metrics_struct.MSSIM = fmt4s(val_mssim);
                    metrics_struct.SAM   = fmt4s(val_sam);
                    metrics_struct.ParamText = param_savetext; 
                    
                    if isempty(summary_metrics)
                        summary_metrics = metrics_struct;
                    else
                        summary_metrics(end+1) = metrics_struct;
                    end
                end
            end
            
            % --- 集計と保存 ---
            if ~isempty(summary_metrics)
                % テーブル変換
                T = struct2table(summary_metrics);
                
                % MPSNRで降順ソート (Bestが上に来るように)
                T = sortrows(T, 'MPSNR', 'descend');
                
                % CSV保存 (全件保存)
                csv_filename = fullfile(dir_result_folder, "summary_metrics.csv");
                writetable(T, csv_filename);
                fprintf("Saved full summary to: %s\n", csv_filename);
                
                % Best Params (.mat) の保存
                best_params_savetext = T.ParamText(1); 
                save(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
                
                % --- 上位 num_show 件の表示 ---
                fprintf("*** Top %d Parameter Sets (by MPSNR) ***\n", num_show);
                disp(T(1:min(num_show, height(T)), :));
                
            else
                warning('No results found in: %s', dir_result_folder);
            end
            
        end
    end
end