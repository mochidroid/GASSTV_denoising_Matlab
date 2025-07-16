clear
addpath(genpath('sub_functions'))

%% Setting parameters
%%%%%%%%%%%%%%%%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blocksize  = [4,4];
% shiftstep  = [1,1];
isTF       = 0;

image = 'Pavia'; % image set = {'Beltsville','Pavia', 'MoffettField'};
deg.sigma       = 0.1; % standard derivation of Gaussian noise
deg.sp_rate     = 0.1; % ratio of salt-and-pepper noise
deg.dl_rate     = 0; % ratio of dead line noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading HSI
image_file_name_tmp = image + "_MixedNoise_" + ...
    "g" + num2str(deg.sigma) + "_sp" + num2str(deg.sp_rate) + "_dl" + num2str(deg.dl_rate) + ...
    "_size_64_64";
image_file_name = "./images/" + image_file_name_tmp + ".mat";
load(image_file_name)

fprintf('~~~ SETTINGS ~~~\n');
fprintf('Image: %s\n', image);
fprintf('sigma: %0.2f s_p: %0.2f\n', deg.sigma, deg.sp_rate);

[n1, n2, n3] = size(org_HSI);

%% Calculating regularization ratio
% SSTV_cir_regularization_ratio = Calc_SSTV_cir_function_value(noisy_HSI) / Calc_SSTV_cir_function_value(org_HSI);
% SSTV_neu_regularization_ratio = Calc_SSTV_neu_function_value(noisy_HSI) / Calc_SSTV_neu_function_value(org_HSI);
% % ASTV_cir_regularization_ratio = Calc_ASTV_cir_function_value(noisy_HSI, blocksize, shiftstep, isTF) ...
% %                                     / Calc_ASTV_cir_function_value(org_HSI, blocksize, shiftstep, isTF);
% STSSTV_cir_regularization_ratio = Calc_STSSTV_cir_function_value(noisy_HSI, blocksize, shiftstep, isTF) ...
%                                     / Calc_STSSTV_cir_function_value(org_HSI, blocksize, shiftstep, isTF);
% 
% fprintf('~~~ R(v)/R(u) ~~~\n');
% fprintf('SSTV circulant: \t%f\n', SSTV_cir_regularization_ratio);
% fprintf('SSTV Neumann:   \t%f\n', SSTV_neu_regularization_ratio);
% % fprintf('ASTV circulant: \t%f\n', ASTV_cir_regularization_ratio);
% fprintf('STSSTV circulant:\t%f\n', STSSTV_cir_regularization_ratio);

%% Calculating STSSTV regularization ratio per blocksize and shiftstep
blocksize_set = [4, 8, 16, 8, 16];
shiftstep_set = [1, 1, 1,  4, 4];

fprintf('~~~ R(v)/R(u) ~~~\n');
for i = 1:numel(blocksize_set)
    blocksize = [blocksize_set(i), blocksize_set(i)];
    shiftstep = [shiftstep_set(i), shiftstep_set(i)];

    tic;
    STSSTV_cir_regularization_ratio = Calc_STSSTV_cir_function_value(noisy_HSI, blocksize, shiftstep, isTF) ...
                                    / Calc_STSSTV_cir_function_value(org_HSI, blocksize, shiftstep, isTF);
    running_time = toc;

    fprintf('blocksize: %d shiftstep: %d ratio: %f time: %f\n', blocksize(1), shiftstep(1), STSSTV_cir_regularization_ratio, running_time);
end
