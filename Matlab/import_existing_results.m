%% import_existing_results.m
% exp_denoising_gstd_Base.m で過去に回した既存法の結果 (dir_save_folder.mat の保存先) を読み込み,
% 古い結果だけを置くフォルダに run_existing_methods.m と同じ階層・同じ形式で保存し直すスクリプト.
%
% 保存先 (run_existing_methods.m の result/ を result/old_results/ に置き換えたもの):
%   <リポジトリ直下>/result/old_results/GASSTV_result/denoising_<image>/<ノイズ条件>/<method>/<params>.mat
%   <リポジトリ直下>/result/old_results/GASSTV_comp_result/...  (復元画像と指標のみ)
%   <リポジトリ直下>/result/old_results/import_log.csv          (全件の判定)
% 比較表は compare_methods.py --existing_dir result/old_results で作れる.
%
% - 対象は run_existing_methods.m が作るはずのファイル (同じ cfg で指定). 保存先に既にあるものは飛ばす.
% - 古い結果は次の 2 点が一致するときだけ取り込む (一致しなければ取り込まずに報告する):
%     1. 観測 HSI_noisy と GT HSI_clean が dataset/<image>/<ノイズ条件>.mat と一致する
%     2. params (epsilon, alpha, beta を含む) が calc_noise_radii.m で計算した値と一致する
% - HSI_clean, HSI_noisy, hsi, deg は dataset のもの, params と復元結果・指標は古い結果のものを保存する.
%
% 実行方法 (リポジトリ直下で):
%   matlab -batch "cfg.noise = 1:9; cfg.rhos = {0.95}; run('Matlab/import_existing_results.m')"
%   cfg の意味は run_existing_methods.m と同じ (noise, images, rhos, methods).
%   cfg.dry_run = true なら保存せずに判定だけ行う.

if ~exist("cfg", "var"), cfg = struct(); end
clearvars -except cfg
close all;

dir_matlab = fileparts(mfilename("fullpath"));
cd(dir_matlab)

addpath(genpath("sub_functions"))
addpath("func_metrics")

fprintf("******* initium *******\n");

%% Conditions (run_existing_methods.m と同じ)
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % 1: g0.1 ps0 pt0
    {0.05,  0.05,  0,     0,    0},     ... % 2: g0.05 ps0.05 pt0 pd0
    {0.1,   0.05,  0,     0,    0},     ... % 3: g0.1 ps0.05 pt0 pd0
    {0.05,  0,     0.05,  0.5,  0},     ... % 4: g0.05 ps0 pt0.05
    {0.1,   0,     0.05,  0.5,  0},     ... % 5: g0.1 ps0 pt0.05
    {0.05,  0.05,  0.05,  0.5,  0},     ... % 6: g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % 7: g0.1 ps0.05 pt0.05 pd0
    {0.1,   0,     0,     0,    0.01},  ... % 8: g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % 9: g0.1 ps0.05 pt0.05 pd0.01
};
idc_noise_conditions = getcfg(cfg, "noise", 1:numel(noise_conditions));

images = {"JasperRidge", "PaviaU"};
idc_images = getcfg(cfg, "images", 1:numel(images));

dir_dataset = fullfile(dir_matlab, "..", "dataset");
dir_result_old = fullfile(dir_matlab, "..", "result", "old_results");  % 古い結果だけを置くフォルダ

load("dir_save_folder.mat", "dir_save_folder");  % 古い結果 (全変数) の保存先

dry_run = getcfg(cfg, "dry_run", false);

rhos = getcfg(cfg, "rhos", {0.90, 0.95});
stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;
maxiter = 20000;


%% Methods (run_existing_methods.m と同じパラメータと保存名)
methods_info = struct([]);

methods_info(end+1).name = "SSTV";
methods_info(end).param_names = {"maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("r%.2f_stop1e-%d", params.rho_radius, stopcri_idx);

HSSTV_omega = {0.01, 0.03, 0.05};
methods_info(end+1).name = "HSSTV_L1";
methods_info(end).param_names = {"L", "omega", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {{"L1"}, HSSTV_omega, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx);

methods_info(end+1).name = "HSSTV_L12";
methods_info(end).param_names = {"L", "omega", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {{"L12"}, HSSTV_omega, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx);

l0l1HTV_stepsize_reduction = {0.999, 0.9999};
l0l1HTV_L10ball_th = {0.02, 0.03};
methods_info(end+1).name = "l0l1HTV";
methods_info(end).param_names = {"L10ball_th", "stepsize_reduction", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {l0l1HTV_L10ball_th, l0l1HTV_stepsize_reduction, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("sr%.5g_th%.2f_r%.2f_maxiter%d", ...
        params.stepsize_reduction, params.L10ball_th, params.rho_radius, maxiter);

LRTDTV_tau = 1;
LRTDTV_lambda_param = {10, 15, 20, 25, sqrt(100*100)};
LRTDTV_lambda = cellfun(@(x) 100 * x / sqrt(100*100), LRTDTV_lambda_param);
LRTDTV_rank = {[100*0.8, 100*0.8, 10], [51, 51, 10]};
methods_info(end+1).name = "LRTDTV";
methods_info(end).param_names = {"tau", "lambda", "rank"};
methods_info(end).params = {LRTDTV_tau, LRTDTV_lambda, LRTDTV_rank};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("l%.4g_r%d_stop1e-%d", params.lambda, params.rank(1), stopcri_idx);

methods_info(end+1).name = "TPTV";
methods_info(end).param_names = {"Rank", "initial_rank", "maxIter", "lambda"};
methods_info(end).params = {{[7,7,5]}, {2}, {50, 100}, {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2}};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("maxiter%d_l%.4g", params.maxIter, params.lambda);

methods_info(end+1).name = "FastHyMix";
methods_info(end).param_names = {"k_subspace"};
methods_info(end).params = {{4, 8, 12}};
methods_info(end).get_params_savetext = @(params) sprintf("sub%d", params.k_subspace);

methods_to_run = getcfg(cfg, "methods", ["SSTV", "HSSTV_L1", "HSSTV_L12", "l0l1HTV"]);
methods_info = methods_info(arrayfun(@(m) any(string(m.name) == string(methods_to_run)), methods_info));


%% Importing
% status: exists (既にある) / imported (取り込んだ) / no_old (古い結果が無い)
%         params_differ (古い結果のパラメータが違う) / obs_differ (観測が違う) / read_error (読めない)
records = cell(0, 6);  % {image, noise, method, setting, status, note}

for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
deg = struct();
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};

name_noise = append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
    "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate));
data = load(fullfile(dir_dataset, image, append(name_noise, ".mat")), "HSI_clean", "HSI_noisy", "hsi", "deg");
HSI_clean = single(data.HSI_clean);
HSI_noisy = single(data.HSI_noisy);
hsi = data.hsi;
deg = data.deg;

for idx_method = 1:numel(methods_info)
name_method = methods_info(idx_method).name;
params_name = methods_info(idx_method).param_names;
[params_comb, num_params_comb] = ParamsList2Comb(methods_info(idx_method).params);

for idx_params_comb = 1:num_params_comb
params = struct();
for idx_params = 1:numel(params_name)
    params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
    if strcmp(params_name{idx_params}, "rho_radius")
        params = calc_noise_radii(params, deg, hsi, HSI_noisy);
    end
end
name_params_savetext = methods_info(idx_method).get_params_savetext(params);

dir_save_method_folder = fullfile(dir_result_old, "GASSTV_result", ...
    append("denoising_", image), name_noise, name_method);
dir_save_comp_method_folder = fullfile(dir_result_old, "GASSTV_comp_result", ...
    append("denoising_", image), name_noise, name_method);
path_save = fullfile(dir_save_method_folder, append(name_params_savetext, ".mat"));
path_old = fullfile(dir_save_folder, append("denoising_", image), name_noise, name_method, ...
    append(name_params_savetext, ".mat"));
row = {image, name_noise, name_method, name_params_savetext};

if isfile(path_save)
    records(end+1, :) = [row, {"exists", ""}]; %#ok<SAGROW>
    continue
end
if ~isfile(path_old)
    records(end+1, :) = [row, {"no_old", ""}]; %#ok<SAGROW>
    continue
end

% H: (Google ドライブ) からの読み込みは失敗することがあるので数回試す
old = [];
for k_try = 1:3
    try
        old = load(path_old);
        break
    catch err
        pause(5);
    end
end
if isempty(old)
    records(end+1, :) = [row, {"read_error", string(err.message)}]; %#ok<SAGROW>
    continue
end

if ~isequal(old.HSI_noisy, HSI_noisy) || ~isequal(old.HSI_clean, HSI_clean)
    records(end+1, :) = [row, {"obs_differ", "HSI_noisy / HSI_clean が dataset と違う"}]; %#ok<SAGROW>
    continue
end
note = compare_params(params, old.params);
if note ~= ""
    records(end+1, :) = [row, {"params_differ", note}]; %#ok<SAGROW>
    continue
end

% 保存 (run_existing_methods.m と同じ変数と形式). params は古い結果のもの (現行の式の値と一致を確認済み)
params             = old.params;
HSI_restored       = old.HSI_restored;
removed_noise      = old.removed_noise;
other_result       = old.other_result;
val_mpsnr          = old.val_mpsnr;
val_mssim          = old.val_mssim;
val_sam            = old.val_sam;
vals_psnr_per_band = old.vals_psnr_per_band;
vals_ssim_per_band = old.vals_ssim_per_band;

% 指標が復元画像から再計算した値と合うかも確認
if abs(calc_MPSNR(HSI_restored, HSI_clean) - val_mpsnr) > 1e-6
    records(end+1, :) = [row, {"metric_differ", "val_mpsnr が再計算値と違う"}]; %#ok<SAGROW>
    continue
end

if ~dry_run
    if ~exist(dir_save_method_folder, "dir"), mkdir(dir_save_method_folder); end
    if ~exist(dir_save_comp_method_folder, "dir"), mkdir(dir_save_comp_method_folder); end
    save(path_save, ...
        "HSI_clean", "HSI_noisy", "hsi", "deg", "image", ...
        "HSI_restored", "removed_noise", ...
        "val_mpsnr", "val_mssim", "val_sam", ...
        "vals_psnr_per_band", "vals_ssim_per_band", ...
        "params", "other_result", ...
        "-v7.3", "-nocompression" ...
    );
    save(fullfile(dir_save_comp_method_folder, append(name_params_savetext, ".mat")), ...
        "HSI_restored", "params", "val_mpsnr", "val_mssim", "val_sam", ...
        "-v7.3", "-nocompression" ...
    );
end
records(end+1, :) = [row, {"imported", sprintf("MPSNR %.4f", val_mpsnr)}]; %#ok<SAGROW>
fprintf("imported: %s\n", path_save);

end
end
end
end


%% Report
T = cell2table(records, "VariableNames", ["image", "noise", "method", "setting", "status", "note"]);
fprintf("\n~~~ SUMMARY ~~~\n");
disp(groupsummary(T, "status"));
missing = T(~ismember(T.status, ["exists", "imported"]), :);
if height(missing) > 0
    fprintf("\n取り込めなかったもの (古い結果が無い / 条件が違う):\n");
    disp(missing);
end
if ~dry_run
    if ~exist(dir_result_old, "dir"), mkdir(dir_result_old); end
    writetable(T, fullfile(dir_result_old, "import_log.csv"), "Encoding", "UTF-8");
end

fprintf("******* finis *******\n");


%% Local functions
function v = getcfg(cfg, name, default)
if isfield(cfg, name)
    v = cfg.(name);
else
    v = default;
end
end


function note = compare_params(p_new, p_old)
% p_new の各フィールドが p_old と一致するか. 違えばその内容を返す ("" なら一致)
note = "";
for fn = string(fieldnames(p_new))'
    if ~isfield(p_old, fn)
        note = note + fn + " が無い; ";
        continue
    end
    a = p_new.(fn); b = p_old.(fn);
    if isnumeric(a) && isnumeric(b)
        same = isequal(size(a), size(b)) && all(abs(double(a(:)) - double(b(:))) <= 1e-6 * max(1, abs(double(a(:)))));  % 古い結果は single で計算したものがある
    else
        same = isequal(string(a), string(b));
    end
    if ~same
        note = note + sprintf("%s: 現行 %s / 古い %s; ", fn, mat2str(double(a), 6), mat2str(double(b), 6));
    end
end
end
