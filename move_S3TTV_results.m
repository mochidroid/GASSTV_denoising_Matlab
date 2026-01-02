%% migrate_old_results_to_new_rule.m
% Old rule:
%   old_root/denoising_<image>/<method>/<name_params_savetext>/image_result.mat
%   old_root/denoising_<image>/<method>/<name_params_savetext>/other_result.mat
%
% New rule (same as exp_denoising_real_Base.m):
%   dir_save_folder/denoising_<image>/<method_new>/<name_params_savetext>.mat
%   dir_save_comp_folder/denoising_<image>/<method_new>/<name_params_savetext>.mat

clear; clc;

%% ====== USER SETTINGS ======
old_root = "H:\マイドライブ\Matlab_sto\S3TTV_for_JSTARS\result";  % ★ここだけ環境に合わせて調整

images = ["IndianPines", "Suwannee"];

methods_old = ["SSTV", "HSSTV_L1", "HSSTV_L12", "l0l1HTV", "STV", "SSST", "FastHyMix"];
methods_new = ["SSTVc", "HSSTV_L1_cir", "HSSTV_L12_cir", "l0l1HTV_cir", "STV_cir", "SSST_cir", "FastHyMix"];


% 現行の保存先（exp_denoising_real_Base.m と同じ前提）
load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

fprintf("=== Migration start ===\n");
fprintf("Old root: %s\n", old_root);
fprintf("New root (main): %s\n", string(dir_save_folder));
fprintf("New root (comp): %s\n", string(dir_save_comp_folder));

%% ====== MAIN LOOP ======
n_ok = 0; n_skip = 0; n_fail = 0;

for iImg = 1:numel(images)
    image = images(iImg);
    old_img_dir = fullfile(old_root, "denoising_" + image);

    if ~isfolder(old_img_dir)
        warning("[SKIP] Image folder not found: %s", old_img_dir);
        n_skip = n_skip + 1;
        continue;
    end

    for iM = 1:numel(methods_old)
        m_old = methods_old(iM);
        m_new = methods_new(iM);

        % 旧 method フォルダ候補（念のため _cir が付いている旧データにも対応）
        cand1 = fullfile(old_img_dir, m_old);
        cand2 = fullfile(old_img_dir, m_new); % もし旧側にも _cir が付いていた場合

        if isfolder(cand1)
            old_method_dir = cand1;
        elseif isfolder(cand2)
            old_method_dir = cand2;
        else
            warning("[SKIP] Method folder not found: %s / %s", cand1, cand2);
            n_skip = n_skip + 1;
            continue;
        end

        % name_params_savetext は「旧 method フォルダ直下のサブフォルダ名」
        d = dir(old_method_dir);
        d = d([d.isdir]);
        subnames = string({d.name});
        subnames = subnames(subnames ~= "." & subnames ~= "..");

        if isempty(subnames)
            warning("[SKIP] No parameter subfolders found in: %s", old_method_dir);
            n_skip = n_skip + 1;
            continue;
        end

        % 新しい保存先 method フォルダ
        new_method_dir_main = fullfile(dir_save_folder, "denoising_" + image, m_new);
        new_method_dir_comp = fullfile(dir_save_comp_folder, "denoising_" + image, m_new);
        if ~isfolder(new_method_dir_main), mkdir(new_method_dir_main); end
        if ~isfolder(new_method_dir_comp), mkdir(new_method_dir_comp); end

        for iP = 1:numel(subnames)
            name_params_savetext = subnames(iP);

            old_param_dir = fullfile(old_method_dir, name_params_savetext);
            f_img   = fullfile(old_param_dir, "image_result.mat");
            f_other = fullfile(old_param_dir, "other_result.mat");

            if ~isfile(f_img) || ~isfile(f_other)
                warning("[SKIP] Missing mat files: %s (img:%d other:%d)", ...
                    old_param_dir, isfile(f_img), isfile(f_other));
                n_skip = n_skip + 1;
                continue;
            end

            try
                Simg = load(f_img);
                Soth = load(f_other);

                % 旧データに依存する前提の変数を取り出し（無ければエラーにする）
                required_img = ["HSI_noisy","hsi","image","HSI_restored","removed_noise"];
                required_oth = ["params","other_result"];

                for rn = required_img
                    if ~isfield(Simg, rn)
                        error("Missing variable '%s' in %s", rn, f_img);
                    end
                end
                for rn = required_oth
                    if ~isfield(Soth, rn)
                        error("Missing variable '%s' in %s", rn, f_other);
                    end
                end

                % まとめて新ルールで保存（現行 script と同じ変数名構成）
                HSI_noisy     = Simg.HSI_noisy;
                hsi           = Simg.hsi;
                image_loaded  = Simg.image; %#ok<NASGU>
                HSI_restored  = Simg.HSI_restored;
                removed_noise = Simg.removed_noise;

                % HSI_clean / hsi 由来の違いがあっても、ここでは新ルールが要求するものだけに寄せる
                params        = Soth.params;
                other_result  = Soth.other_result;

                out_main = fullfile(new_method_dir_main, name_params_savetext + ".mat");
                out_comp = fullfile(new_method_dir_comp, name_params_savetext + ".mat");

                save(out_main, ...
                    "HSI_noisy","hsi","image","HSI_restored","removed_noise","params","other_result", ...
                    "-v7.3","-nocompression");

                save(out_comp, ...
                    "HSI_restored","params", ...
                    "-v7.3","-nocompression");

                n_ok = n_ok + 1;
                fprintf("[OK] %s | %s -> %s\n", image, m_new, name_params_savetext);

            catch ME
                n_fail = n_fail + 1;
                warning("[FAIL] %s | %s | %s : %s", image, m_new, name_params_savetext, ME.message);
            end
        end
    end
end

fprintf("\n=== Migration summary ===\n");
fprintf("OK   : %d\n", n_ok);
fprintf("SKIP : %d\n", n_skip);
fprintf("FAIL : %d\n", n_fail);
fprintf("=== Migration end ===\n");
