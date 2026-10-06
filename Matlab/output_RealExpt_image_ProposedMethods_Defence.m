clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

%% Switching operator
% is_show_HSI = 1;  
is_output_image = 1; % 画像保存フラグ

%% Setting common parameters
images = {...
    "IndianPines", ...
    "Suwannee", ...
};

% Radius parameters; 1st: IndianPines, 2nd: Suwannee
epsilon = {30, 100};
alpha = {200, 800};
beta = {100, 5000};

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;

% 保存先フォルダパスのロード
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

% --- 追加設定: Classification用 ---
train_ratio = 0.10; % 学習データの割合 (10%)
rnd_seed = 50;      % 再現性確保のためのシード
% Ground Truthのパス
gt_path = 'H:\マイドライブ\MATLAB_Share\HSIData\IndianPine\Indian_pines_gt.mat'; 
% ---------------------------------


%% Setting output parameters
idx_output = 1;

switch idx_output
    case 1
        idx_image = 1; % IndianPines

        gain_restoration = 1.5; 
        gain_diff = 5;

        save_band = 32; % IndianPines

        crop_start_pos = [2, 30];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';

        arrow_head_pos = [20, 8];
        % arrow_head_pos = [38, 74];
        % arrow_head_pos = [53, 64];
        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_dir_tblr = "l";
        arrow_methods_idc = [];
        
    case 2
        idx_image = 2; % Suwannee

        gain_restoration = 1;
        gain_diff = 8;
        save_band = 196; 
        
        crop_start_pos = [30, 80];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'bl';
end

 
%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "output_name", "SSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% l0l1HTV
l0l1HTV.stepsize_reduction = 0.999; 
l0l1HTV.L10ball_th = 0.03; 

methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "output_name", "$\llHTV$", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_sr%.5g_th%.2f_maxiter%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            l0l1HTV.stepsize_reduction, l0l1HTV.L10ball_th, maxiter), ...
    "line_style", "--", ...
    "enable", false ...
);

% HSSTV_L1
HSSTV.omega = 0.05;

methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "output_name", "HSSTV1", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "output_name", "HSSTV2", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% TPTV
TPTV.Rank = [7,7,5];
TPTV.initial_rank = 2;
TPTV.maxIter = 50; 
TPTV.lambdas = 1e-3;

methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "output_name", "TPTV", ...
    "get_params_savetext", ...
        sprintf("maxiter%d_l%.4g", TPTV.maxIter, TPTV.lambdas), ...
    "line_style", "--", ...
    "enable", true ...
);

% QRNN3D
QRNN3D.ckpt = "complex"; 

methods_info(end+1) = struct( ...
    "name", "QRNN3D", ...
    "output_name", "QRNN3D", ...
    "get_params_savetext", ...
        sprintf("p_%s_0normalized", QRNN3D.ckpt), ...
    "line_style", "--", ...
    "enable", true ...
);


% GASSTV_CondatVu
GASSTV_CondatVu.sigma_sp = "90"; 
GASSTV_CondatVu.sigma_l = "med";
GASSTV_CondatVu.lambda_rho_sp = 0.9; 
GASSTV_CondatVu.lambda2 = 0.1; 
GASSTV_CondatVu.num_segments = 4;
GASSTV_CondatVu.order_filt = 5; 

methods_info(end+1) = struct( ...
    "name", "GASSTV_CondatVu", ...
    "output_name", "GASSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_l%.2g_%.2g_sig%s_%s_ns%d_fil%d_stop1e-%d", ...
        epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
        GASSTV_CondatVu.lambda_rho_sp, GASSTV_CondatVu.lambda2, ...
        GASSTV_CondatVu.sigma_sp, GASSTV_CondatVu.sigma_l, ...
        GASSTV_CondatVu.num_segments, GASSTV_CondatVu.order_filt, ...
        stopcri_idx), ...
    "line_style", "-", ...
    "enable", true ...
);

% S3TTV
S3TTV.blocksize = [10, 10];
S3TTV.epsilon = 15;
S3TTV.alpha = 200;
S3TTV.beta = 100;

methods_info(end+1) = struct( ...
    "name", "S3TVT_cir", ...
    "output_name", "S3TTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_bl%d_stop1e-%d", ...
        S3TTV.epsilon, S3TTV.alpha, S3TTV.beta, ...
        S3TTV.blocksize(1), stopcri_idx), ...
    "line_style", "-", ...
    "enable", true ...
);

% GeoSSTV
GeoSSTV.lambda = 0.01;
methods_info(end+1) = struct( ...
    "name", "GeoSSTV", ...
    "output_name", "GeoSSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            20, alpha{idx_image}, beta{idx_image}, ...
            GeoSSTV.lambda, stopcri_idx), ...
    "line_style", "-", ...
    "enable", true ...
);



methods_info = methods_info([methods_info.enable]); % removing true methods
num_methods = numel(methods_info);


%% Loading clean and noisy HS images
image = images{idx_image};

[HSI_noisy, hsi] = Load_real_HSI(image);
HSI_noisy = single(HSI_noisy);

fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);


%% Cropping clean and noisy HS images
% image_noisy = HSI_noisy(:,:,save_band)*gain_restoration;

image_noisy = HSI_noisy(:,:,save_band)*gain_restoration;
image_noisy = Crop_Embed_image(image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

if ~isempty(arrow_methods_idc)
    image_noisy = Embed_arrow(image_noisy, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
end


%% Loading and cropping result images
for idx_method = 1:num_methods 
    name_method = methods_info(idx_method).name;

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        name_method ...
    );

    params_savetext = methods_info(idx_method).get_params_savetext;

    load(fullfile(dir_result_folder, append(params_savetext, ".mat")), ...
        "HSI_restored");

    methods_info(idx_method).HSI_restored = HSI_restored;

    % Cropping result images
    % image_restored = HSI_restored(:,:,save_band)*gain_restoration;
    % 
    % methods_info(idx_method).image_restored         = image_restored;

    image_restored = HSI_restored(:,:,save_band)*gain_restoration;
    image_restored = Crop_Embed_image(image_restored, ...
                crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

    
    if find(arrow_methods_idc == idx_method)
        image_restored = Embed_arrow(image_restored, ...
                    arrow_head_pos, arrow_length, ...
                    arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    end

    methods_info(idx_method).image_restored         = image_restored;


end


%% Showing cropped images
% cat_restored_image = image_noisy;
% 
% for idx_method = 1:num_methods 
%     cat_restored_image = cat(2, cat_restored_image, methods_info(idx_method).image_restored);
% end
% 
% figure;
% imshow(cat_restored_image);
% 
% figure;
% imshow(cat_crop_restored_image);

cat_restored_image = image_noisy;

for idx_method = 1:num_methods 
    cat_restored_image = cat(2, cat_restored_image, methods_info(idx_method).image_restored);
end

figure;
imshow(cat_restored_image);


%% Showing HSI
if exist("is_show_HSI", "var") && is_show_HSI == 1
    cat_HSI = HSI_noisy;

    for idx_method = 1:num_methods
        cat_HSI = cat(2, cat_HSI, methods_info(idx_method).HSI_restored);
    end

    cat_diff_GeoSSTV = abs(repmat(methods_info(num_methods).HSI_restored, [1, num_methods+1, 1]) - cat_HSI) * gain_diff;

    implay(cat(1, cat_HSI, cat_diff_GeoSSTV));
end





%% =========================================================================
%% =========================================================================
%% Output / Evaluation Section: Images & Classification
%% =========================================================================
% NOTE:
% - Classification itself is ALWAYS executed (independent of output flag).
% - Saving restored images / classification maps is controlled by is_output_image.

do_output = (exist("is_output_image", "var") && is_output_image == 1);

% --- Output folder configuration (only when saving) ---
if do_output
    dir_output_root = fullfile(...
        dir_save_comp_folder, ...
        "Defence", ...
        sprintf("%s_b%d", image, save_band));

    dir_output_result_folder      = fullfile(dir_output_root, "restored_image");
    dir_output_class_folder       = fullfile(dir_output_root, "classification"); % classification maps + metrics

    if ~exist(dir_output_result_folder, 'dir'),      mkdir(dir_output_result_folder); end
    if ~exist(dir_output_class_folder, 'dir'),       mkdir(dir_output_class_folder); end
else
    dir_output_root = "";
    dir_output_result_folder = "";
    dir_output_class_folder = "";
end

%% --- Classification preparation (ALWAYS) ---
fprintf('\n--- Preparing Classification ---\n');

% Load ground-truth labels (ALWAYS)
% - If you already define gt_path above, it will be used here.
% - The GT variable name inside the .mat file may differ; we load the first field.
if exist("gt_path", "var") && isfile(gt_path)
    S_gt = load(gt_path);
    f_gt = fieldnames(S_gt);
    GT_original = S_gt.(f_gt{1});
else
    error("GT file not found. Please set gt_path correctly: %s", gt_path);
end

% Crop GT to match the evaluation region (if needed)
if strcmp(image, "IndianPines")
    GT_cropped = GT_original(1:end-25, 26:end);
else
    GT_cropped = GT_original;
end

% Make sure the noisy/restored HSI used for classification has the same spatial size as GT_cropped
if strcmp(image, "IndianPines") && size(HSI_noisy, 1) ~= size(GT_cropped, 1)
    HSI_noisy_cls = HSI_noisy(1:end-25, 26:end, :);
else
    HSI_noisy_cls = HSI_noisy;
end

% Colormaps
cmap_base  = jet(16);
cmap_class = [0 0 0; cmap_base]; % background: black
cmap_error = [0 0 0; 1 0 0];     % black: correct, red: error

% Save Ground Truth map (only when saving)
if do_output
    imwrite(uint8(GT_cropped), cmap_class, fullfile(dir_output_class_folder, "GroundTruth_Map.png"));
end

% Table to store metrics (always)
metrics_table = table();

%% --- 1) Noisy image: classification (ALWAYS) ---
fprintf('Processing Noisy Image ...\n');
train_ratio = 0.1; % <- change if needed

[acc_noisy, map_noisy, err_noisy] = run_svm_classification(HSI_noisy_cls, GT_cropped, train_ratio, rnd_seed);
metrics_table = [metrics_table; table("Noisy", acc_noisy.OA, acc_noisy.AA, acc_noisy.Kappa, ...
    'VariableNames', {'Method', 'OA', 'AA', 'Kappa'})];

fprintf('  -> OA: %.2f%%\n', acc_noisy.OA);

% Save images/maps (only when saving)
if do_output
    % A) Save restored image (2D) of the selected band (noisy)
    imwrite(image_noisy, fullfile(dir_output_result_folder, "image_noisy.png"), "BitDepth", 8);
    % B) Save classification maps (noisy)
    save_classification_maps(map_noisy, err_noisy, "Noisy", dir_output_class_folder, cmap_class, cmap_error);
end

%% --- 2) Restored images: classification (ALWAYS) ---
for idx_method = 1:num_methods
    name_method = methods_info(idx_method).name;
    name_output = methods_info(idx_method).output_name;

    fprintf('Processing %s ...\n', name_method);

    % A) Save restored image(s) (only when saving)
    if do_output
        image_restored = methods_info(idx_method).image_restored;
        imwrite(image_restored, ...
            fullfile(dir_output_result_folder, sprintf("image_%s.png", name_output)), ...
            "BitDepth", 8);
    end

    % B) Run classification (ALWAYS)
    HSI_restored = methods_info(idx_method).HSI_restored;

    if strcmp(image, "IndianPines") && size(HSI_restored, 1) ~= size(GT_cropped, 1)
        HSI_restored_cls = HSI_restored(1:end-25, 26:end, :);
    else
        HSI_restored_cls = HSI_restored;
    end

    [acc, map_class, map_err] = run_svm_classification(HSI_restored_cls, GT_cropped, train_ratio, rnd_seed);

    % Store metrics (always)
    metrics_table = [metrics_table; table(string(name_method), acc.OA, acc.AA, acc.Kappa, ...
        'VariableNames', {'Method', 'OA', 'AA', 'Kappa'})];

    fprintf('  -> OA: %.2f%%\n', acc.OA);

    % Save classification maps (only when saving)
    if do_output
        save_classification_maps(map_class, map_err, name_output, dir_output_class_folder, cmap_class, cmap_error);
    end
end

%% --- Summary (ALWAYS) ---
fprintf('\n=== Classification Summary ===\n');
disp(metrics_table);

% Save summary table (only when saving)
if do_output
    writetable(metrics_table, fullfile(dir_output_class_folder, "metrics_summary.csv"));
    fprintf("\nAll results saved to:\n%s\n", dir_output_root);
end
function [acc, class_map, error_map] = run_svm_classification(HSI, GT, train_ratio, seed)
    rng(seed);
    [nr, nc, nb] = size(HSI);
    X_vector = reshape(HSI, nr*nc, nb);
    gt_vector = double(GT(:));
    
    idx_labeled = find(gt_vector > 0);
    X_labeled = X_vector(idx_labeled, :);
    y_labeled = gt_vector(idx_labeled);
    
    % Z-score Normalization (Important for SVM)
    X_labeled = zscore(X_labeled);
    
    % Split
    cv = cvpartition(y_labeled, 'HoldOut', 1 - train_ratio);
    idx_train = cv.training;
    idx_test = cv.test;
    
    % Train SVM (Linear kernel)
    t = templateSVM('Standardize', true, 'KernelFunction', 'linear');
    Mdl = fitcecoc(X_labeled(idx_train, :), y_labeled(idx_train), ...
        'Learners', t, 'Coding', 'onevsone', 'Options', statset('UseParallel', false));
    
    % Predict
    pred_test = predict(Mdl, X_labeled(idx_test, :));
    
    % Metrics
    conf_mat = confusionmat(y_labeled(idx_test), pred_test);
    oa = trace(conf_mat) / sum(conf_mat(:));
    class_acc = diag(conf_mat) ./ sum(conf_mat, 2);
    aa = mean(class_acc(~isnan(class_acc)));
    po = oa;
    pe = (sum(conf_mat, 1) * sum(conf_mat, 2)) / (sum(conf_mat(:))^2);
    kappa = (po - pe) / (1 - pe);
    
    acc.OA = oa * 100;
    acc.AA = aa * 100;
    acc.Kappa = kappa * 100;
    
    % Create Maps
    pred_all = predict(Mdl, X_labeled);
    class_map = zeros(nr, nc);
    class_map(idx_labeled) = pred_all;
    
    error_map = zeros(nr, nc);
    is_error = (GT > 0) & (GT ~= class_map);
    error_map(is_error) = 1;
end

function save_classification_maps(class_map, error_map, method_name, save_dir, cmap_c, cmap_e)
    f_class = fullfile(save_dir, sprintf("ClassMap_%s.png", method_name));
    imwrite(uint8(class_map), cmap_c, f_class);
    
    f_error = fullfile(save_dir, sprintf("ErrorMap_%s.png", method_name));
    imwrite(uint8(error_map), cmap_e, f_error);
end
