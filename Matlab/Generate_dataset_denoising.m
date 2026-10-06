%% Generate_dataset_denoising.m
% Load_HSI で読み込んだ GT にノイズを重畳し, GT と観測のペアを dataset フォルダに保存する.
% 観測の生成手順は exp_denoising_gstd_Base.m と同じ (noise_seed = "default").
%
% 保存先: <リポジトリ直下>/dataset/<image>/g<..>_ps<..>_pt<..>_pd<..>.mat
%   HSI_clean : GT (single, [0,1] 正規化済み)
%   HSI_noisy : 観測 (single)
%   hsi       : 画像サイズなどの情報
%   deg       : ノイズ条件と各ノイズ成分 (gaussian_noise, sparse_noise, stripe_noise, deadline_mpose)
%   image     : 画像名

clear
close all;

dir_matlab = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(dir_matlab, "sub_functions")))

fprintf("******* initium *******\n");

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

idc_noise_conditions = 1:numel(noise_conditions);

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

idc_images = 1:numel(images);

dir_dataset = fullfile(dir_matlab, "..", "dataset");

% true にすると既存のファイルも作り直す
overwrite = false;


%% Generating dataset
for idx_image = idc_images
for idx_noise_condition = idc_noise_conditions
    deg = struct();
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
    image = images{idx_image};

    name_noise = append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate));
    dir_save = fullfile(dir_dataset, image);
    path_save = fullfile(dir_save, append(name_noise, ".mat"));

    if ~overwrite && isfile(path_save)
        fprintf("skip (exists): %s\n", path_save);
        continue
    end

    % exp_denoising_gstd_Base.m と同じ手順
    [HSI_clean, hsi] = Load_HSI(image);
    noise_seed = "default";
    [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

    HSI_clean = single(HSI_clean);
    HSI_noisy = single(HSI_noisy);

    % ノイズ成分も single で保存 (容量削減)
    for fn = ["gaussian_noise", "sparse_noise", "stripe_noise", "deadline_mpose"]
        if isfield(deg, fn)
            deg.(fn) = single(deg.(fn));
        end
    end

    if ~exist(dir_save, "dir"), mkdir(dir_save); end
    save(path_save, "HSI_clean", "HSI_noisy", "hsi", "deg", "image", "-v7");
    fprintf("saved: %s  (size: %d x %d x %d)\n", path_save, hsi.n1, hsi.n2, hsi.n3);
end
end

fprintf("******* finis *******\n");
