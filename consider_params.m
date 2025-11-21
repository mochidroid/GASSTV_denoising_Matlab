clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")

%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     tint  pd
    {0.1,   0,     0,     0,    0},     ... % g0.1 ps0 pt0
    {0.1,   0.05,  0,     0,    0},     ... % g0.1 ps0.05 pt0 pd0
    {0.1,   0,     0.05,  0.5,  0},     ... % g0.1 ps0 pt0.05
    {0.1,   0,     0,     0,    0.01},  ... % g0.1 ps0 pt0 pd0.01
    {0.1,   0.05,  0.05,  0.5,  0.01},  ... % g0.1 ps0.05 pt0.05 pd0.01
    {0.05,  0.05,  0.05,  0.5,  0},     ... % g0.05 ps0.05 pt0.05 pd0
    {0.1,   0.05,  0.05,  0.5,  0},     ... % g0.1 ps0.05 pt0.05 pd0
};

idx_noise_condition = 2;

images = {...
    "JasperRidge", ...
    "PaviaU", ...
    "Beltsville", ...
};

idx_image = 2;

fmt4s = @(x) round(x,4,"significant");

name_method = 'GASSTV_GLR_Const_OraGuide';
% name_method = 'GASSTV_Oraguide_Const';
% name_method = 'GASSTV_GLR_Oraguide';



%% Setting condition
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};


load("dir_save_comp_folder.mat", "dir_save_comp_folder");

dir_result_folder = fullfile(...
    dir_save_comp_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);


%% 
load(fullfile(dir_result_folder, "summary_metrics.mat"));


% "param_names", {{"lambda_rho_sp", "lambda2", "sigma_sp", ...
%         "num_segments", "maxiter", "stopcri", "rho_radius"}}, ...

% label1 = "lambda_rho_sp";
% label2 = "lambda_rho_l";

label1 = "lambda_rho_sp";
label2 = "num_segments";


%% 
%% summary_metrics から値を取り出し
% m1  = [summary_metrics.lambda1];
% m2 = [summary_metrics.lambda2];

m1  = [summary_metrics.(label1)];
m2 = [summary_metrics.(label2)];

% m1  = [summary_metrics.sigma_sp];
% m2 = [summary_metrics.sigma_s];
% m2 = [summary_metrics.sigma_l];

zz  = [summary_metrics.mpsnr];

%% ユニークな組み合わせを求めて平均を取る
[groups, ~, idx] = unique([m1(:), m2(:)], 'rows');
avg_mpsnr = accumarray(idx, zz(:), [], @mean);

ux = unique(m1);
uy = unique(m2);
Z  = nan(numel(ux), numel(uy));
for k = 1:size(groups,1)
    i = find(ux == groups(k,1));
    j = find(uy == groups(k,2));
    Z(i,j) = avg_mpsnr(k);
end

%% 3D surface
figure;
surf(ux, uy, Z');
shading interp; grid on; colorbar;
xlabel(label1);
ylabel(label2);
zlabel('MPSNR (dB, averaged)');
title('Average MPSNR vs \lambda_{\rho,sp}, \lambda_{\rho,l}');

%% 2D heatmap も欲しければ
figure;
imagesc(ux, uy, Z');
set(gca,'YDir','normal'); colorbar;
xlabel(label1); ylabel(label2);
title('Average MPSNR (dB)');