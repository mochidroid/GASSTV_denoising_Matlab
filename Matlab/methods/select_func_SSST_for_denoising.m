function [HSI_restored, removed_noise, other_result] ...
     = select_func_SSST_for_denoising(HSI_clean, HSI_noisy, params, deg)

addpath('./methods/SSST');

% Selecting SSST function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_SSST_g_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_SSST_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_SSST_gt_for_denoising(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_SSST_gst_for_denoising(HSI_clean, HSI_noisy, params);
end