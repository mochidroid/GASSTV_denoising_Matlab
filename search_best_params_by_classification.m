clear;
close all;
addpath(genpath("sub_functions"))
addpath("func_metrics")
fprintf("******* Real Data Classification Benchmark *******\n");

%% 1. Config & Data Selection
% -------------------------------------------------------------------------
image_name = "IndianPines"; % 評価したい画像
% image_name = "Suwannee"; 

% Indian PinesのGround Truthパス (環境に合わせて変更してください)
gt_path = 'H:\マイドライブ\MATLAB_Share\HSIData\IndianPine\Indian_pines_gt.mat';

% 結果フォルダのルート
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

% 分類設定
train_ratio = 0.10; % 10%学習, 90%テスト
rnd_seed = 42;      % 再現性確保のため固定

%% 2. Loading Ground Truth
% -------------------------------------------------------------------------
if isfile(gt_path)
    load(gt_path, 'indian_pines_gt');
    GT_original = indian_pines_gt;
else
    error('GT file not found. Please check the path.');
end

% Indian Pines特有のクロップ処理 (復元画像とサイズを合わせる)
if strcmp(image_name, "IndianPines")
    GT_cropped = GT_original(1:end-25, 26:end);
else
    GT_cropped = GT_original; % 他の画像ならそのまま(必要に応じて調整)
end
fprintf("GT Loaded. Size: %dx%d\n", size(GT_cropped,1), size(GT_cropped,2));

%% 3. Setting Method Parameters (GASSTV_CondatVu)
% -------------------------------------------------------------------------
% 比較したいパラメータ範囲を設定
GASSTV_CondatVu.sigma_sp = {"90", "med"};
GASSTV_CondatVu.sigma_l = {"med"};
GASSTV_CondatVu.lambda_rho_sp = [0.9, 1, 1.1];
GASSTV_CondatVu.lambda2 = [0.1, 1, 10];
GASSTV_CondatVu.num_segments = [4];
GASSTV_CondatVu.order_filt = [5, 1];

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;
maxiter = 20000;

epsilon = 30;
alpha = 200;
beta = 100;

% メソッド情報構造体
methods_info(1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", "sigma_l", ...
        "num_segments", "order_filt", "maxiter", "stopcri", ...
        "epsilon", "alpha", "beta"}}, ...
    "params", {{GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, maxiter, stopcri, epsilon, alpha, beta}}, ...
    "get_params_savetext", @(params) ...
        sprintf("e%g_a%g_b%g_l%.2g_%.2g_sig%s_%s_ns%d_fil%d_stop1e-%d", ...
        params.epsilon, params.alpha, params.beta, ...
        params.lambda_rho_sp, params.lambda2, params.sigma_sp, params.sigma_l, ...
        params.num_segments, params.order_filt, stopcri_idx), ...
    "enable", true ...
);

methods_info = methods_info([methods_info.enable]); 
num_methods = numel(methods_info);

%% 4. Main Loop: Load Denoised Image -> Classify -> Record
% -------------------------------------------------------------------------
for idx_method = 1:num_methods
    name_method = methods_info(idx_method).name;
    params_name = methods_info(idx_method).param_names;
    params_cell = methods_info(idx_method).params;
    [params_comb, num_params_comb] = ParamsList2Comb(params_cell);

    % 保存先フォルダ (Simと違いRealはフォルダ階層が浅いことが多いので注意)
    % 既存コードのフォルダ構成: denoising_IndianPines/MethodName/
    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image_name), ...
        name_method ...
    );
    
    fprintf("\n~~~ METHOD: %s (%d combinations) ~~~\n", name_method, num_params_comb);
    fprintf("Folder: %s\n", dir_result_folder);
    
    % 結果格納用
    summary_metrics = struct([]);
    
    % パラメータループ
    for idx_params_comb = 1:num_params_comb
        % 1. パラメータの取得とファイル名生成
        params = struct();
        for idx_params = 1:numel(params_name)
            params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
        end
        
        param_savetext = methods_info(idx_method).get_params_savetext(params);
        mat_filename = fullfile(dir_result_folder, append(param_savetext, ".mat"));
        
        % 2. 復元画像のロード
        if isfile(mat_filename)
            load(mat_filename, "HSI_restored");
            
            % 3. SVM分類の実行 (関数化)
            fprintf("[%d/%d] Classifying: %s ... ", idx_params_comb, num_params_comb, param_savetext);
            [OA, AA, Kappa] = run_classification(HSI_restored, GT_cropped, train_ratio, rnd_seed);
            fprintf("OA: %.4f, AA: %.4f, Kappa: %.4f\n", OA, AA, Kappa);
            
            % 4. 結果の保存
            % パラメータ情報をコピー
            metrics_struct = params; 
            % 評価指標を追加
            metrics_struct.OA = OA;
            metrics_struct.AA = AA;
            metrics_struct.Kappa = Kappa;
            
            if isempty(summary_metrics)
                summary_metrics = metrics_struct;
            else
                summary_metrics(end+1) = metrics_struct;
            end
        else
            fprintf("[%d/%d] File not found: %s\n", idx_params_comb, num_params_comb, param_savetext);
        end
    end
    
    % 5. テーブル化と保存
    if ~isempty(summary_metrics)
        T = struct2table(summary_metrics);
        
        % OAでソート (降順)
        T = sortrows(T, 'OA', 'descend');
        
        % CSV保存
        csv_filename = fullfile(dir_result_folder, "summary_classification_metrics.csv");
        writetable(T, csv_filename);
        fprintf("\nResult table saved to: %s\n", csv_filename);
        
        % 最良結果の表示
        fprintf("\n*** Best Parameter Set (by OA) ***\n");
        disp(T(1, :));
        
        % 全表示
        disp(T);
    else
        warning("No results were processed. Check file paths.");
    end
end


%% =========================================================================
%% Helper Function: Run SVM Classification
%% =========================================================================
function [OA, AA, Kappa] = run_classification(HSI, GT, train_ratio, seed)
    % 再現性確保
    rng(seed);
    
    % データ整形
    [nr, nc, nb] = size(HSI);
    X_vector = reshape(HSI, nr*nc, nb);
    gt_vector = double(GT(:));
    
    % ラベルありデータのみ抽出
    idx_labeled = find(gt_vector > 0);
    X_labeled = X_vector(idx_labeled, :);
    y_labeled = gt_vector(idx_labeled);
    
    % 正規化 (Z-score) - SVMの精度に重要
    X_labeled = zscore(X_labeled);
    
    % 学習・テスト分割
    cv = cvpartition(y_labeled, 'HoldOut', 1 - train_ratio);
    idx_train = cv.training;
    idx_test = cv.test;
    
    X_train = X_labeled(idx_train, :);
    y_train = y_labeled(idx_train);
    X_test = X_labeled(idx_test, :);
    y_test = y_labeled(idx_test);
    
    % SVM学習 (Linear Kernel for speed and stability)
    t = templateSVM('Standardize', true, 'KernelFunction', 'linear');
    % 'Quiet', true で学習中のログを抑制
    Mdl = fitcecoc(X_train, y_train, 'Learners', t, 'Coding', 'onevsone', ...
        'Options', statset('UseParallel', false)); 
    
    % 予測
    pred_test = predict(Mdl, X_test);
    
    % 指標計算
    conf_mat = confusionmat(y_test, pred_test);
    
    % OA
    OA = trace(conf_mat) / sum(conf_mat(:));
    
    % AA
    class_acc = diag(conf_mat) ./ sum(conf_mat, 2);
    valid_mask = ~isnan(class_acc);
    AA = mean(class_acc(valid_mask));
    
    % Kappa
    po = OA;
    pe = (sum(conf_mat, 1) * sum(conf_mat, 2)) / (sum(conf_mat(:))^2);
    Kappa = (po - pe) / (1 - pe);
end