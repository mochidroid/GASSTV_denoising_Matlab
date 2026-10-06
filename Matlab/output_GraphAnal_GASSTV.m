%% output_graph_validity_10images_v3.m
% - guideimage_*.png : imwrite (no margin)
% - spagraph_*.png   : imwrite (no margin)
% - spegraph1..4.png : white background + white legend

clear; close all; clc;
addpath(genpath("./sub_functions"));

%% ===== Selecting conditions =====
noise_conditions = { ...
    {0.1,   0,     0,     0,    0},     ... % case1
    {0.05,  0.05,  0,     0,    0},     ...
    {0.1,   0.05,  0,     0,    0},     ...
    {0.05,  0,     0.05,  0.5,  0},     ...
    {0.1,   0,     0.05,  0.5,  0},     ...
    {0.05,  0.05,  0.05,  0.5,  0},     ...
    {0.1,   0.05,  0.05,  0.5,  0},     ... % case7
};

idx_case1 = 1;
idx_case7 = 7;

images = {"JasperRidge","PaviaU","Beltsville"};
idx_image = 1;
image = images{idx_image};

%% ===== Parameters =====
sigma_sp_opts = "90";
w_idx = 2;

num_segments = 4;
k_lap        = 10;
sigma_l      = "med";
order_filt   = 5;

% dpi_png = 600;
dpi_png = 300;

is_output = 1;

%% ===== Load save directory (YOUR SPECIFICATION) =====
load("dir_save_comp_folder.mat", "dir_save_comp_folder");

save_figure_folder = fullfile( ...
    dir_save_comp_folder, ...
    "GASSTV_Defence", ...
    "ParamAnal");

if exist("is_output","var")
    if ~exist(save_figure_folder,"dir")
        mkdir(save_figure_folder);
    end
end

%% ===== Load data =====
[U, hsi] = Load_HSI(image);
U = single(U);

labels = {"clean","case1","case7"};
V_cases = { ...
    U, ...
    gen_noisy(U, noise_conditions{idx_case1}), ...
    gen_noisy(U, noise_conditions{idx_case7}) ...
};

%% ============================================================
%  (1) guideimage_*.png  (imwrite, no margin)
%% ============================================================
for t = 1:3
    guide = mean(V_cases{t}, 3);
    guide = guide - min(guide(:));
    guide = guide / max(guide(:));
    guide_u8 = uint8(255 * guide);

    if exist("is_output","var")
        fname = fullfile(save_figure_folder, ...
            sprintf("guideimage_%s.png", labels{t}));
        imwrite(guide_u8, fname, "BitDepth", 8);
        fprintf("[Saved] %s\n", fname);
    end
end

%% ============================================================
%  (2) spagraph_*.png  (imwrite, no margin)
%% ============================================================
for t = 1:3
    [Wsp, ~] = Create_SpatialGraphWeight(V_cases{t}, sigma_sp_opts);
    Wsp_slice = squeeze(Wsp(:,:,1,:));
    W = Wsp_slice(:,:,w_idx);

    cmap = parula(256);
    W_u8 = uint8(255 * min(max(W,0),1));
    W_rgb = ind2rgb(W_u8+1, cmap);

    if exist("is_output","var")
        fname = fullfile(save_figure_folder, ...
            sprintf("spagraph_%s.png", labels{t}));
        imwrite(W_rgb, fname);
        fprintf("[Saved] %s\n", fname);
    end
end

%% ============================================================
%  (3) spegraph1..4.png  (white bg + WHITE legend)
%% ============================================================
R_all = cell(1,3);
for t = 1:3
    [~, ~, ~, info] = Create_SpectralGraphLaplacian( ...
        V_cases{t}, num_segments, sigma_l, k_lap, order_filt);
    R_all{t} = info.representative;
end

Splot = min(num_segments, size(R_all{1},1));

col_clean = [0 0 0];
col_case1 = [0.85 0 0];
col_case7 = [0 0.35 1.0];

for s = 1:Splot
    fig = figure("Visible","off", "Position",[100 100 1100 650], ...
                 "Color","w");

    plot(R_all{1}(s,:), '-', 'Color', col_clean, 'LineWidth', 3.0); hold on;
    plot(R_all{2}(s,:), '-', 'Color', col_case1, 'LineWidth', 3.0);
    plot(R_all{3}(s,:), '-', 'Color', col_case7, 'LineWidth', 3.0);
    hold off;

    grid on; axis tight;
    set(gca, ...
        "Color","w", ...
        "XColor","k", "YColor","k", ...
        "LineWidth",1.4, "FontSize",20);

    xlabel("Band index", "Color","k");
    ylabel("Intensity", "Color","k");

    %%% CHANGED: white legend
    lgd = legend({"clean","case1","case7"}, "Location","best");
    set(lgd, ...
        "TextColor","k", ...
        "Color","w", ...
        "EdgeColor","k");

    % title(sprintf("Representative spectrum (segment %d)", s), "Color","k");

    if exist("is_output","var")
        fname = fullfile(save_figure_folder, ...
            sprintf("spegraph%d.png", s));
        print(fig, fname, "-dpng", sprintf("-r%d", dpi_png));
        fprintf("[Saved] %s\n", fname);
    end
    close(fig);
end

fprintf("******* finis *******\n");

%% ===== helper =====
function V = gen_noisy(U, cond)
    deg.gaussian_sigma      = cond{1};
    deg.sparse_rate         = cond{2};
    deg.stripe_rate         = cond{3};
    deg.stripe_intensity    = cond{4};
    deg.deadline_rate       = cond{5};
    noise_seed = "default";
    [V, ~] = Generate_obsv_for_denoising(U, deg, noise_seed);
    V = single(V);
end
