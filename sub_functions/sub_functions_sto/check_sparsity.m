clear
close all;

addpath('function_S3TTV')
addpath(genpath('sub_functions'))

rng('default')

%% Generating observation
%%%%%%%%%%%%%%%%%%%%% User settings of experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg.is_gaussian         = 1;
deg.is_sparse           = 1;
deg.is_stripe           = 0;

deg.gaussian_sigma      = 0.1; % standard derivation of Gaussian noise

deg.sparse_rate         = 0.05;

deg.stripe_rate         = 0.05;
deg.stripe_intensity    = 0.5;


image = 'JasperRidge';
% image = 'PaviaU';

[HSI_clean, hsi] = Load_HSI(image);
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

%% 
U_per = permute(HSI_clean, [2,3,1]);
V_per = permute(HSI_noisy, [2,3,1]);

DbU = Db_Neumann(HSI_clean);
DbV = Db_Neumann(HSI_noisy);

DbDvV = Dv_Neumann(DbV);
DbDvU = Dv_Neumann(DbU);

DbV_per = permute(DbV, [2,3,1]);
DbU_per = permute(DbU, [2,3,1]);

DbDvV_per = permute(DbDvV, [2,3,1]);
DbDvU_per = permute(DbDvU, [2,3,1]);


%% 
implay(cat(1, cat(2, U_per, V_per), (cat(2, DbU_per, DbV_per)+0.5)./2, (cat(2, DbDvU_per, DbDvV_per)+0.5)./2))
