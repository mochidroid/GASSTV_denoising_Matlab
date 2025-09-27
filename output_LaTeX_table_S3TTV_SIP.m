clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.1 pt0 pd0
};

idc_noise_conditions = 1;
num_conditions = numel(idc_noise_conditions);


images = {...
    "JasperRidge64", ...
    "PaviaU64", ...
};

images_outputs = {...
    'Jasper Ridge', ...
    'Pavia University', ...
};

idc_images = 1:numel(images);
% idc_images = [2, 1];
num_images = numel(idc_images);


%% Setting common parameters
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

dir_save_LaTeX_txt = fullfile(dir_save_comp_folder, "S3TTVSPDHG_SIP");
mkdir(dir_save_LaTeX_txt)

 
%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "output_name", "SSTV", ...
    "enable", true ...
);

% l0l1HTV
methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "output_name", "$\llHTV$", ...
    "enable", true ...
);

% HSSTV_L1
methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "output_name", "HSSTV1", ...
    "enable", true ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "output_name", "HSSTV2", ...
    "enable", true ...
);


% GeoSSTV
methods_info(end+1) = struct( ...
    "name", "S3TTV_NPPDS", ...
    "output_name", "S3TTV NPPDS", ...
    "enable", true ...
);

% GeoSSTV
methods_info(end+1) = struct( ...
    "name", "S3TTV_SPDHG", ...
    "output_name", "S3TTV SPDHG", ...
    "enable", true ...
);



methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Loading HS images and Calculating metrics
vals_mpsnr = nan(num_methods, num_conditions, num_images);
vals_mssim = nan(num_methods, num_conditions, num_images);

for idx_image = idc_images
for idx_noise_condition = idc_noise_conditions
    %% Selecting noise condition
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
    image = images{idx_image};
    
    
    [HSI_clean, hsi] = Load_HSI(image);
    noise_seed = "default";
    [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);
    
    HSI_clean = single(HSI_clean);
    HSI_noisy = single(HSI_noisy);
    
    
    %% Loading results
    for idx_method = 1:num_methods
        name_method = methods_info(idx_method).name;
    
        dir_result_folder = fullfile(...
            dir_save_comp_folder, ...
            append("denoising_", image), ...
            append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                    "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
            name_method ...
        );
    
        load(fullfile(dir_result_folder, "best_params.mat"), "best_params_savetext");
        methods_info(idx_method).get_params_savetext = best_params_savetext;
    
        load(fullfile(dir_result_folder, append(best_params_savetext, ".mat")), ...
            "val_mpsnr", "val_mssim");
    
        vals_mpsnr(idx_method, idx_noise_condition, idx_image) = val_mpsnr;
        vals_mssim(idx_method, idx_noise_condition, idx_image) = val_mssim;
    end

end
end


%% Outputting LaTeX Table
num_cols   = 3 + num_methods;
cmid_all   = sprintf('\\cmidrule(lr){1-%d}', num_cols);
cmid_inner = sprintf('\\cmidrule(lr){2-%d}', num_cols);

% 数値の表示桁
fmt_mpsnr = '%.2f';
fmt_mssim = '%.4f';

% 本文バッファ（char）
text_tex = '';

out_txt = fullfile(dir_save_LaTeX_txt, 'table_MPSNR_MSSIM_rows.txt');
fid = fopen(out_txt,'w');

for idx_image = idc_images
    image_output = images_outputs{idx_image};

    % ここからは全て char の連結 [] と sprintf で改行を入れる
    text_tex = [text_tex, cmid_all, sprintf(' \n')];
    text_tex = [text_tex, sprintf('\\multirow{%d}{*}{%s} \n', 2*num_conditions, image_output)];

    for idx_noise_condition = idc_noise_conditions
        % 現ケースの値ベクトル（methods_info の順番のまま）
        mpsnr_vec = squeeze(vals_mpsnr(:, idx_noise_condition, idx_image)).';
        mssim_vec = squeeze(vals_mssim(:, idx_noise_condition, idx_image)).';

        % ===== MPSNR 行：最良(太字)・次点(下線)を判定 =====
        mpsnr_cells = format_with_rank_char(mpsnr_vec, fmt_mpsnr); % cell array of char
        mpsnr_join  = strjoin(mpsnr_cells, ' & ');                  % char

        % 最初のケースだけは "Case" を同じ行で、以降は内側の cmidrule を入れる
        if idx_noise_condition == idc_noise_conditions(1)
            text_tex = [text_tex, sprintf('& \\multirow{2}{*}{Case %d} & \n', idx_noise_condition)];
        else
            text_tex = [text_tex, cmid_inner, sprintf(' \n')];
            text_tex = [text_tex, sprintf('& \\multirow{2}{*}{Case %d} & \n', idx_noise_condition)];
        end
        text_tex = [text_tex, 'MPSNR & ', mpsnr_join, sprintf(' \\\\ \n')];

        % ===== MSSIM 行：最良(太字)・次点(下線)を判定 =====
        mssim_cells = format_with_rank_char(mssim_vec, fmt_mssim);
        mssim_join  = strjoin(mssim_cells, ' & ');
        text_tex = [text_tex, '& & MSSIM & ', mssim_join, sprintf(' \\\\ \n')];
    end

    text_tex = [text_tex, sprintf('\n')];
    if idx_image == idc_images(end)
        text_tex = [text_tex, sprintf('\\bottomrule\n')];
    end
end

% ファイル出力（text_tex は既に改行文字を含む char）
fprintf(fid, '%s', text_tex);
fclose(fid);
fprintf('LaTeX rows written to: %s\n', out_txt);


%% ===== 補助関数（char 版）=====
function cell_strs = format_with_rank_char(vals, fmt)
% vals: 1×M 数値ベクトル（NaN 可）
% fmt : sprintf 用のフォーマット（例 '%.2f'）
% 役割: 最大値を \textbf{}、2番目（同率含む）を \underline{} で装飾し char の cell 配列で返す

    M = numel(vals);
    cell_strs = cell(1, M);

    if all(isnan(vals))
        [cell_strs{:}] = deal('--');
        return;
    end

    finite_vals = vals(~isnan(vals));
    vmax = max(finite_vals);
    second_candidates = finite_vals(finite_vals < vmax);
    has_second = ~isempty(second_candidates);
    if has_second
        vsecond = max(second_candidates);
    else
        vsecond = NaN;
    end

    eps_tie = 1e-12;
    for k = 1:M
        if isnan(vals(k))
            cell_strs{k} = '--';
            continue;
        end
        numtxt = sprintf(fmt, vals(k));
        if abs(vals(k) - vmax) <= eps_tie
            cell_strs{k} = ['\textbf{', numtxt, '}'];
        elseif has_second && abs(vals(k) - vsecond) <= eps_tie
            cell_strs{k} = ['\underline{', numtxt, '}'];
        else
            cell_strs{k} = numtxt;
        end
    end
end
