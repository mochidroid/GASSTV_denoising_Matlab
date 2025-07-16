function [HSI_restored, removed_noise, other_result] ...
     = select_func_RISSTV_for_denoising(HSI_clean, HSI_noisy, params, deg)

% Selecting RISSTV function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_RISSTV_g_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_RISSTV_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_RISSTV_gt_for_denoising(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_RISSTV_gst_for_denoising(HSI_clean, HSI_noisy, params);
end