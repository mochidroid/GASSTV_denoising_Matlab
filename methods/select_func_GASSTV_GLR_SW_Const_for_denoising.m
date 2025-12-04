function [HSI_restored, removed_noise, other_result] ...
     = select_func_GASSTV_GLR_SW_Const_for_denoising(HSI_clean, HSI_noisy, params, deg)

addpath("./methods/GASSTV_GLR_SW/")

% Selecting GASSTV_OraGuide function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_GLR_SW_Const_g_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_GLR_SW_Const_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_GLR_SW_Const_gt_for_denoising(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_GASSTV_GLR_SW_Const_gst_for_denoising(HSI_clean, HSI_noisy, params);
end