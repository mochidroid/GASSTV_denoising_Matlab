function [HSI_restored, removed_noise, other_result] ...
     = select_func_GASSTV_for_denoising(HSI_clean, HSI_noisy, params, deg)

addpath("./methods/GASSTV/")

% Selecting GASSTV function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_g_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_gt_for_denoising(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_gst_for_denoising(HSI_clean, HSI_noisy, params);
end