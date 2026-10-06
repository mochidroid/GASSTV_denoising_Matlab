function [HSI_restored, removed_noise, other_result] ...
     = select_func_l0l1HTV_for_denoising(HSI_clean, HSI_noisy, params, deg)

addpath('./methods/l0l1HTV');

% Selecting l0l1HTV function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_l0l1HTV_g_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_l0l1HTV_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_l0l1HTV_gt_for_denoising(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_l0l1HTV_gst_for_denoising(HSI_clean, HSI_noisy, params);
end