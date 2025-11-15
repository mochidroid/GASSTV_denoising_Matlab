function [HSI_restored, removed_noise, other_result] = func_TPTV(HSI_clean, HSI_noisy, params, ~)
addpath(genpath('./methods/TPTV/'));
param = params;

Ori_H = HSI_clean;
Noi_H = HSI_noisy;

sizeof = size(HSI_noisy);


% Saving current directory
dir_prev = pwd;

cd('./methods/TPTV/')
addpath(genpath('functions'))

[ output_image,U_x,V_x,E,Xo] = WETV(Noi_H,Ori_H, param);

HSI_restored = reshape(output_image, sizeof);

removed_noise.all_noise = reshape(E, sizeof);

other_result.U_x = U_x;
other_result.V_x = V_x;
other_result.Xo = Xo;
    

% Returning current directory
cd(dir_prev)

end

