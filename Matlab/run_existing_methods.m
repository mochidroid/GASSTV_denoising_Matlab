%% run_existing_methods.m
% 既存法 (比較手法) のデノイジング実験を回すスクリプト.
% パラメータは exp_denoising_gstd_Base.m と同じ. 提案法 (GASSTV) は Python 版 (run_gasstv.py) で回す.
%
% - 観測は <リポジトリ直下>/dataset/<image>/<ノイズ条件>.mat を使う
%   (無ければ Generate_dataset_denoising.m と同じ手順でその場で生成する).
% - 結果は Python 版と同じ階層に保存する:
%     <リポジトリ直下>/result/GASSTV_result/denoising_<image>/g<..>_ps<..>_pt<..>_pd<..>/<method>/<params>.mat
%     <リポジトリ直下>/result/GASSTV_comp_result/...  (復元画像と指標のみ)
% - SSTV / HSSTV / l0l1HTV には g, gt 用の関数が無いので, gs, gst 用の関数で代用する.
%   スパースノイズも deadline ノイズも無いときは alpha = 0 となり S = 0 に固定されるため, g, gt 用と等価.
% - epsilon, alpha, beta は calc_noise_radii.m で計算する (deadline ノイズはスパースの項に含める).
% - 既存法のコード (methods/ 以下) には手を加えていない.
%
% 実行方法 (リポジトリ直下で):
%   matlab -batch "run('Matlab/run_existing_methods.m')"
% 設定は下の既定値か, 実行前に cfg を置いて上書きする (指定しない項目は既定値):
%   matlab -batch "cfg.noise = [1 9]; cfg.rhos = {0.90}; cfg.methods = {'SSTV','HSSTV_L1'}; run('Matlab/run_existing_methods.m')"
%   cfg.noise   : ノイズ条件の番号 (下の noise_conditions の行番号)
%   cfg.images  : 画像の番号 (下の images の番号)
%   cfg.rhos    : rho_radius の候補 (cell)
%   cfg.methods : 実行する手法名 (cell または string 配列). SSTV, HSSTV_L1, HSSTV_L12, l0l1HTV, LRTDTV, TPTV, FastHyMix
%   cfg.skip_existing : true なら結果ファイルがある組み合わせを飛ばす

if ~exist("cfg", "var"), cfg = struct(); end
clearvars -except cfg
close all;

dir_matlab = fileparts(mfilename("fullpath"));
cd(dir_matlab)  % FastHyMix などが Matlab/ からの相対パスを使うため

addpath(genpath("sub_functions"))
addpath("func_metrics")
addpath("methods")

fprintf("******* initium *******\n");

%% Selecting conditions
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

idc_noise_conditions = getcfg(cfg, "noise", 9);

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

idc_images = getcfg(cfg, "images", 1:numel(images));

dir_dataset = fullfile(dir_matlab, "..", "dataset");
dir_result  = fullfile(dir_matlab, "..", "result");

% true にすると結果ファイルが既にある組み合わせを飛ばす
skip_existing = getcfg(cfg, "skip_existing", true);


%% Setting common parameters
% rho は提案法 (Python) と共通の値で比べる (0.90 と 0.95 のどちらに固定するかを決めるため両方回す)
rhos = getcfg(cfg, "rhos", {0.90, 0.95});

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;


%% Setting each methods info
% "paths": その手法の実行中だけ追加するフォルダ (genpath で追加)
methods_info = struct([]);

% SSTV
methods_info(end+1).name = "SSTV";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    run_tv_method("SSTV", HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("r%.2f_stop1e-%d", params.rho_radius, stopcri_idx);
methods_info(end).paths = {"methods/SSTV"};
methods_info(end).enable = true;

% HSSTV_L1
HSSTV_omega = {0.01, 0.03, 0.05};
methods_info(end+1).name = "HSSTV_L1";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    run_tv_method("HSSTV", HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"L", "omega", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {{"L1"}, HSSTV_omega, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx);
methods_info(end).paths = {"methods/HSSTV"};
methods_info(end).enable = true;

% HSSTV_L12
methods_info(end+1).name = "HSSTV_L12";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    run_tv_method("HSSTV", HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"L", "omega", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {{"L12"}, HSSTV_omega, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("o%.2f_r%.2f_stop1e-%d", params.omega, params.rho_radius, stopcri_idx);
methods_info(end).paths = {"methods/HSSTV"};
methods_info(end).enable = true;

% l0l1HTV
l0l1HTV_stepsize_reduction = {0.999, 0.9999};
l0l1HTV_L10ball_th = {0.02, 0.03};
methods_info(end+1).name = "l0l1HTV";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    run_tv_method("l0l1HTV", HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"L10ball_th", "stepsize_reduction", "maxiter", "stopcri", "rho_radius"};
methods_info(end).params = {l0l1HTV_L10ball_th, l0l1HTV_stepsize_reduction, maxiter, stopcri, rhos};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("sr%.5g_th%.2f_r%.2f_maxiter%d", ...
        params.stepsize_reduction, params.L10ball_th, params.rho_radius, maxiter);
methods_info(end).paths = {"methods/l0l1HTV"};
methods_info(end).enable = true;

% LRTDTV
LRTDTV_tau = 1;
LRTDTV_lambda_param = {10, 15, 20, 25, sqrt(100*100)};
LRTDTV_lambda = cellfun(@(x) 100 * x / sqrt(100*100), LRTDTV_lambda_param);
LRTDTV_rank = {[100*0.8, 100*0.8, 10], [51, 51, 10]};
methods_info(end+1).name = "LRTDTV";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    func_LRTDTV(HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"tau", "lambda", "rank"};
methods_info(end).params = {LRTDTV_tau, LRTDTV_lambda, LRTDTV_rank};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("l%.4g_r%d_stop1e-%d", params.lambda, params.rank(1), stopcri_idx);
methods_info(end).paths = {"methods/LRTDTV"};
methods_info(end).enable = true;

% TPTV
TPTV_Rank = {[7,7,5]};
TPTV_initial_rank = {2};
TPTV_maxIter = {50, 100};
TPTV_lambdas = {5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2};
methods_info(end+1).name = "TPTV";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    func_TPTV_for_denoising(HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"Rank", "initial_rank", "maxIter", "lambda"};
methods_info(end).params = {TPTV_Rank, TPTV_initial_rank, TPTV_maxIter, TPTV_lambdas};
methods_info(end).get_params_savetext = @(params) ...
    sprintf("maxiter%d_l%.4g", params.maxIter, params.lambda);
methods_info(end).paths = {"methods/TPTV"};
methods_info(end).enable = true;

% FastHyMix (func_FastHyMix 内で methods/HSI-MixedNoiseRemoval-FastHyMix-main に cd する)
% single で渡すと内部で NaN になるため double で渡す (run_fasthymix)
% なお JasperRidge のスパース+ストライプ条件では FastHyMix 本体 (.p) が内部エラーを出し,
% func_FastHyMix が "Error" を表示して観測をそのまま返す (FastHyMix 側の問題)
FastHyMix_k_subspace = {4, 8, 12};
methods_info(end+1).name = "FastHyMix";
methods_info(end).func = @(HSI_clean, HSI_noisy, params, deg) ...
    run_fasthymix(HSI_clean, HSI_noisy, params, deg);
methods_info(end).param_names = {"k_subspace"};
methods_info(end).params = {FastHyMix_k_subspace};
methods_info(end).get_params_savetext = @(params) sprintf("sub%d", params.k_subspace);
methods_info(end).paths = {};
methods_info(end).enable = true;


% 実行する手法 (既定は rho を使う 4 手法. LRTDTV, TPTV, FastHyMix は rho に依存しないので必要なときだけ)
methods_to_run = getcfg(cfg, "methods", ["SSTV", "HSSTV_L1", "HSSTV_L12", "l0l1HTV"]);
for k = 1:numel(methods_info)
    methods_info(k).enable = any(string(methods_info(k).name) == string(methods_to_run));
end

methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Running Expt.
idx_exp = 0;
total_exp = numel(idc_noise_conditions) * numel(idc_images);

for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
%% Loading observation
deg = struct();
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};

name_noise = append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
    "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate));
path_dataset = fullfile(dir_dataset, image, append(name_noise, ".mat"));

if isfile(path_dataset)
    data = load(path_dataset, "HSI_clean", "HSI_noisy", "hsi", "deg");
    HSI_clean = data.HSI_clean;
    HSI_noisy = data.HSI_noisy;
    hsi = data.hsi;
    deg = data.deg;
    fprintf("Observation: %s\n", path_dataset);
else
    % Generate_dataset_denoising.m と同じ手順で生成
    [HSI_clean, hsi] = Load_HSI(image);
    [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, "default");
    fprintf("Observation: generated (%s が無いため)\n", path_dataset);
end

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

idx_exp = idx_exp + 1;


%% Running methods
for idx_method = 1:num_methods
name_method = methods_info(idx_method).name;
func_method = methods_info(idx_method).func;
params_name = methods_info(idx_method).param_names;
params_cell = methods_info(idx_method).params;

[params_comb, num_params_comb] = ParamsList2Comb(params_cell);

% 手法ごとのパスを追加
method_paths = cellfun(@(p) genpath(fullfile(dir_matlab, p)), ...
    methods_info(idx_method).paths, "UniformOutput", false);
for k = 1:numel(method_paths), addpath(method_paths{k}); end

for idx_params_comb = 1:num_params_comb

params = struct();
for idx_params = 1:numel(params_name)
    params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};

    % If rho_radius exists, calculate epsilon, alpha, and beta (deadline はスパースの項に含める)
    if strcmp(params_name{idx_params}, "rho_radius")
        params = calc_noise_radii(params, deg, hsi, HSI_noisy);
    end
end

name_params_savetext = methods_info(idx_method).get_params_savetext(params);

dir_save_method_folder = fullfile(dir_result, "GASSTV_result", ...
    append("denoising_", image), name_noise, name_method);
dir_save_comp_method_folder = fullfile(dir_result, "GASSTV_comp_result", ...
    append("denoising_", image), name_noise, name_method);
path_save = fullfile(dir_save_method_folder, append(name_params_savetext, ".mat"));

fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Method: %s\n", name_method);
fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Noise: %s (stripe intensity %g)\n", name_noise, deg.stripe_intensity);
fprintf("Parameter settings: %s\n", name_params_savetext)
fprintf("Methods: (%d/%d), Cases: (%d/%d), Params:(%d/%d)\n", ...
    idx_method, num_methods, idx_exp, total_exp, idx_params_comb, num_params_comb);

if skip_existing && isfile(path_save)
    fprintf("-> skip (exists)\n");
    continue
end

rng("default");
[HSI_restored, removed_noise, other_result] ...
    = func_method(HSI_clean, HSI_noisy, params, deg);


% Calculating metrics
val_mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
val_mssim  = calc_MSSIM(HSI_restored, HSI_clean);
val_sam    = calc_SAM(HSI_restored, HSI_clean);

fprintf("~~~ RESULTS ~~~\n");
fprintf("MPSNR: %#.4g\n", val_mpsnr);
fprintf("MSSIM: %#.4g\n", val_mssim);
fprintf("SAM  : %#.4g\n", val_sam);

[vals_psnr_per_band, vals_ssim_per_band] = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


% Saving each result
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
fprintf("Saved: %s\n", path_save);

close all

end

for k = 1:numel(method_paths), rmpath(method_paths{k}); end

end
end
end

fprintf("******* finis *******\n");


%% Local functions
function v = getcfg(cfg, name, default)
% cfg に name があればその値, なければ default
if isfield(cfg, name)
    v = cfg.(name);
else
    v = default;
end
end


function [HSI_restored, removed_noise, other_result] = run_tv_method(name, HSI_clean, HSI_noisy, params, deg)
% SSTV / HSSTV / l0l1HTV の呼び出し.
% g, gt 用の関数は存在しないため, ストライプ無しなら gs 用, 有りなら gst 用を使う
% (スパースノイズも deadline ノイズも無いときは alpha = 0 で S = 0 に固定されるので等価).
if deg.stripe_rate == 0
    func = str2func(append("func_", name, "_gs_for_denoising"));
else
    func = str2func(append("func_", name, "_gst_for_denoising"));
end
[HSI_restored, removed_noise, other_result] = func(HSI_clean, HSI_noisy, params);
end


function [HSI_restored, removed_noise, other_result] = run_fasthymix(HSI_clean, HSI_noisy, params, deg)
% FastHyMix は single 入力だと NaN を出す (func_FastHyMix 内の try/catch で "Error" となり
% 観測がそのまま返る) ため, double に変換して呼び, 出力を single に戻す.
[HSI_restored, removed_noise, other_result] = ...
    func_FastHyMix(double(HSI_clean), double(HSI_noisy), params, deg);
HSI_restored = single(HSI_restored);
if isstruct(removed_noise)
    removed_noise.all_noise = single(removed_noise.all_noise);
end
end
